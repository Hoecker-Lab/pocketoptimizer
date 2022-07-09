import itertools
import multiprocessing as mp
import os
from functools import partial
from typing import List, Dict, Union, Tuple, NoReturn
import logging
import numpy as np
from ffevaluation.ffevaluate import FFEvaluate
from moleculekit.molecule import Molecule
from tqdm.auto import tqdm

logging.root.handlers = []
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - [%(levelname)s] - %(message)s",
    handlers=[
        logging.FileHandler(os.environ.get('POCKETOPTIMIZER_LOGFILE')),
        logging.StreamHandler()
    ]
)

logger = logging.getLogger('pocketoptimizer.scoring.ffscorer')


class FFScorer:
    """
    FFScorer class for ligand protein interactions (ligand sidechain and ligand scaffold)
    Utilizing molecular dynamics binding energy related computations
    Only taking into account non-bonded interactions (VdW and electrostatic terms)
    Requires a built protein ligand complex, as well as complexes for each single mutation.
    Each ligand pose is scored against fixed scaffold and against each sidechain rotamer respectively

    """
    def __init__(self, ligand: str, ligand_poses: str, project_path: str, forcefield: str, rotamer_path: str):
        """
        Initialize a :class: FFScorer object

        Parameters
        ----------
        ligand: str
            Path to ligand mol2 file
        ligand_poses: str
            Path to ligand poses file
        project_path: str
            Path to project/working directory
        forcefield: str
            Force field used for scoring
        rotamer_path: str
            Path to sampled side chain rotamers
        """
        self.ligand = ligand
        self.ligand_poses = ligand_poses
        self.project_path = project_path
        self.forcefield = forcefield
        self.rotamer_path = rotamer_path
        self.lig_scaff = ''
        self.lig_side = ''

    @staticmethod
    def calculate_ff(poses: Tuple[int], struc: Molecule, ffev: FFEvaluate, pose_coords: np.ndarray = None,
                     rot_coords: np.ndarray = None, resid: str = None, chain: str = None) -> Tuple[np.float]:
        """Calculates the force field components and returns the energies

        Parameters
        ----------
        pose: tuple
            Tuple of ligand pose and rotamer id
        struc: :class: moleculekit.Molecule
            Object of Molecule to set the rotamer coordinates
        ffev: :class:ffevaluation.ffevaluate.FFevaluate
            Object to calculate the energies from between sets ligand and side chain/scaffold
        rot_coords: list
            Rotamer coordinates np.array(natoms, 3, nrots)
        pose_coords: list
              Ligand pose coordinates np.array(natoms, 3, nposes)
        resid: str
            Residue id of the rotamer

        Returns
        -------
        Tuple containing vdW and electrostatic energy for the scored ligand pose
        """
        if len(poses) == 2:
            struc.set('coords', rot_coords[..., poses[1]], f'chain {chain} and resid {resid}')

        # set coordinates to current ligand pose
        struc.set('coords', pose_coords[..., poses[0]], 'segid L')
        energies = ffev.calculateEnergies(struc.coords)
        nrg = energies['vdw'], energies['elec']
        return nrg

    def run_self_nrg(self, pose_coords: np.ndarray, mutations: List[Dict[str, Union[str, List[str]]]], nposes: int, outpath: str, ncpus: int = 1) -> NoReturn:
        """
        Calculates self energies

        Parameters
        ----------
        pose_coords: np.ndarray
            ligand pose coordinates in a numpy array
        mutations: List[Dict]
            List of dictionaries defining the mutations at different positions
        nposes: int
            Number of conformers
        outpath: str
            Output directory
        ncpus: int
            number of CPUs to use for scoring [default: 1]
        """
        from pocketoptimizer.utility.utils import load_ff_parameters, write_energies, calculate_chunks

        structure_path = os.path.join(self.project_path, 'scaffold', self.forcefield, 'protein_params', 'native_complex')

        if not os.path.exists(os.path.join(structure_path, 'structure.pdb')):
            logger.error('Missing protein-ligand complex build. Run prepare_protein() and include ligand in binding pocket.')
            raise FileNotFoundError('Missing protein-ligand complex build. Run prepare_protein() and include the ligand in the binding pocket.')

        outfile = os.path.join(outpath, 'ligand.csv')

        # We only want to exclude sidechain atoms
        selection = 'protein and not ('

        for position in mutations:
            selection += f"chain {position['chain']} and resid {position['resid']} and sidechain or "

        selection = selection[0:-4] + ')'
        struc, prm = load_ff_parameters(structure_path=structure_path, forcefield=self.forcefield)
        ffev = FFEvaluate(struc, prm, betweensets=(selection, 'segid L'))
        self_nrgs = np.zeros((nposes, 2))

        with mp.Pool(processes=ncpus) as pool:
            with tqdm(total=nposes, desc='Ligand/Scaffold') as pbar:
                for index, energy in enumerate(pool.imap(
                        partial(
                            self.calculate_ff,
                            struc=struc.copy(),
                            ffev=ffev,
                            pose_coords=pose_coords
                        ), [(pose,) for pose in range(nposes)],
                        chunksize=calculate_chunks(nposes=nposes, ncpus=ncpus))):
                    self_nrgs[index] = energy
                    pbar.update()

        # Save data as csv
        write_energies(outpath=outfile,
                       energies=self_nrgs,
                       energy_terms=['VdW', 'ES'],
                       name_a='ligand_pose',
                       nconfs_a=nposes)

    def run_pair_nrg(self, pose_coords: np.ndarray, mutations: List[Dict[str, Union[str, List[str]]]], nposes: int, ncpus: int = 1) -> NoReturn:
        """
        Calculate sidechain ligand energies

        Parameters
        ----------
        pose_coords: np.ndarray
            Ligand pose coordinates
        mutations: list
            List of dictionaries defining the mutations at different positions
        nposes: int
            Number of ligand poses
        ncpus: int
            Number of Cpu's to use for scoring [default: 1]
        """
        from pocketoptimizer.utility.utils import load_ff_parameters, write_energies, calculate_chunks

        structure_path_prepend = os.path.join(self.project_path, 'scaffold', self.forcefield, 'protein_params')
        output_file_prepend = os.path.join(self.project_path, 'energies', f'{self.forcefield}_{self.rotamer_path.split("/")[-1]}', f'ligand_sidechain_{self.forcefield}')
        for position in mutations:
            chain, resid = position['chain'], position['resid']
            rotamer_path_prepend = os.path.join(self.rotamer_path, f'{chain}_{resid}')
            for resname in position['mutations']:
                output_file = os.path.join(output_file_prepend, f'ligand_{chain}_{resid}_{resname}.csv')
                # Check if pairwise-interaction energy file already exists
                if os.path.isfile(output_file):
                    logger.info(f'Ligand-Sidechain-Interaction-Energies with Residue: {chain}_{resid}_{resname} already computed.')
                    continue
                logger.info(
                    f'Ligand-Sidechain-Interaction-Energies with Residue: {chain}_{resid}_{resname} not computed yet.')
                # Create Molecule object from rotamers of each residue
                # Check if rotamers are computed for all mutations
                try:
                    rotamer_mol = Molecule(os.path.join(rotamer_path_prepend, f'{resname}.pdb'))
                except:
                    logger.error(f'Missing rotamer for residue: {chain}_{resid}_{resname}.')
                    raise FileNotFoundError(
                        f'Missing rotamer for residue: {chain}_{resid}_{resname}.')
                nrots = rotamer_mol.coords.shape[-1]
                structure_path = os.path.join(structure_path_prepend, f'{chain}_{resid}_{resname}')

                if not os.path.exists(os.path.join(structure_path, 'structure.pdb')):
                    logger.error(f'Missing mutated protein-ligand-complex build for mutation {chain}_{resid}_{resname}.')
                    raise FileNotFoundError(f'Missing mutated protein-ligand-complex build for mutation {chain}_{resid}_{resname}.')

                struc, prm = load_ff_parameters(structure_path=structure_path, forcefield=self.forcefield)
                ffev = FFEvaluate(struc, prm,
                                  betweensets=(f'chain {chain} and resid {resid} and sidechain', 'segid L'))

                batch = list(itertools.product(*[np.arange(nposes), np.arange(nrots)]))
                # Create array (nposes * nrots * 2 (vdw and es))
                pair_nrgs = np.zeros((nposes, nrots * 2))

                with tqdm(total=len(batch), desc=f'Ligand/{chain}_{resid}_{resname}') as pbar:
                    with mp.Pool(processes=ncpus) as pool:
                        for index, energy in enumerate(pool.imap(
                                partial(
                                    self.calculate_ff,
                                    struc=struc.copy(),
                                    ffev=ffev,
                                    pose_coords=pose_coords,
                                    rot_coords=rotamer_mol.coords,
                                    chain=chain,
                                    resid=resid
                                ), batch,
                                chunksize=calculate_chunks(nposes=nposes, ncpus=ncpus))):
                            pair_nrgs[batch[index][0], batch[index][1]*2:batch[index][1]*2 + 2] = energy
                            pbar.update()

                # Save data as csv
                write_energies(outpath=output_file,
                               energies=pair_nrgs,
                               energy_terms=['VdW', 'ES'],
                               name_a='ligand_pose',
                               nconfs_a=nposes,
                               name_b=resname,
                               nconfs_b=nrots)

    def run_ff_scorer(self, mutations: List[Dict[str, Union[str, List[str]]]], ncpus: int = 1) -> NoReturn:
        """
        Run Ligand ff scoring

        Parameters
        ----------
        mutations: list
            List of dictionaries defining the mutations at different positions
        ncpus: int
            Number of CPUs to use for scoring [default: Number of CPUs available]
        """
        logger.info(f'Score ligand interactions using {self.forcefield} force field.')

        # Output ligand scaffold energy
        self.lig_scaff = os.path.join(self.project_path, 'energies', f'{self.forcefield}_{self.rotamer_path.split("/")[-1]}', f'ligand_scaffold_{self.forcefield}')
        os.makedirs(self.lig_scaff, exist_ok=True)
        # Output ligand side chain energies
        self.lig_side = os.path.join(self.project_path, 'energies', f'{self.forcefield}_{self.rotamer_path.split("/")[-1]}', f'ligand_sidechain_{self.forcefield}')
        os.makedirs(self.lig_side, exist_ok=True)

        try:
            ligand_poses = Molecule(self.ligand_poses)
            ligand_poses.read(f'{os.path.splitext(self.ligand_poses)[0]}.xtc')
        except:
            logger.error('Missing ligand poses.')
            raise FileNotFoundError('Missing ligand poses.')

        pose_coords = ligand_poses.coords
        nposes = pose_coords.shape[-1]

        if os.path.isfile(os.path.join(self.lig_scaff, 'ligand.csv')):
            logger.info(f'Ligand-Scaffold-Interaction-Energies already computed.')

        else:
            logger.info(f'Ligand-Scaffold-Interaction-Energies not computed yet.')
            self.run_self_nrg(pose_coords=pose_coords, mutations=mutations, ncpus=ncpus, nposes=nposes, outpath=self.lig_scaff)

        self.run_pair_nrg(pose_coords=pose_coords, mutations=mutations, nposes=nposes, ncpus=ncpus)
        logger.info('Ligand scoring was successful.')
