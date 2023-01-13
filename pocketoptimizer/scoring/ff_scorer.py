import itertools
import multiprocessing as mp
import os
from functools import partial
from typing import Tuple, NoReturn
import logging
import numpy as np
from ffevaluation.ffevaluate import FFEvaluate
from moleculekit.molecule import Molecule
from tqdm.auto import tqdm

from pocketoptimizer.utility.utils import Storer

logger = logging.getLogger(__name__)


class LigandScorer(Storer):
    """
    Ligand scoring class for ligand protein interactions (ligand sidechain and ligand scaffold)
    Utilizing molecular dynamics binding energy related computations
    Requires a built protein ligand complex, as well as complexes for each single mutation.

    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.energy_terms = ['Vdw', 'Elec']

    def calculate_energy(self, ids: Tuple[int], struc: Molecule, ffev: FFEvaluate, pose_coords: np.ndarray,
                        rot_coords: np.ndarray = None, chain: str = None, resid: str = None) -> np.ndarray:
        """Calculates the force field components and returns the energies

        Parameters
        ----------
        ids: tuple
            Tuple of ligand and rotamer id
        struc: :class: moleculekit.Molecule
            Object of Molecule to set the rotamer coordinates
        ffev: :class:ffevaluation.ffevaluate.FFevaluate
            Object to calculate the energies from between sets ligand and side chain/scaffold
        pose_coords: np.ndarray
              Ligand pose coordinates np.array(natoms, 3, nposes)
        rot_coords: np.ndarray
            Rotamer coordinates np.array(natoms, 3, nrots)
        chain: str
            Chain Id of rotamer
        resid: str
            Residue Id of rotamer

        Returns
        -------
        Array containing vdW and electrostatic energy for the scored ligand pose
        """
        if len(ids) == 2:
            struc.set('coords', rot_coords[..., ids[1]], f'chain {chain} and resid {resid}')

        # set coordinates to current ligand pose
        struc.set('coords', pose_coords[..., ids[0]], 'segid L')
        energies = ffev.calculateEnergies(struc.coords)
        return np.array([energies['vdw'], energies['elec']])

    def run_self_nrg(self, pose_coords: np.ndarray, nposes: int) -> NoReturn:
        """
        Calculates self energies

        Parameters
        ----------
        pose_coords: np.ndarray
            ligand pose coordinates in a numpy array
        nposes: int
            Number of conformers
        """
        from pocketoptimizer.utility.utils import load_ff_parameters, write_energies, calculate_chunks
        from pocketoptimizer.utility.molecule_types import backbone_atoms

        structure_path = os.path.join(self.work_dir, 'scaffold', self.forcefield, 'protein_params', 'ligand_sampling_pocket')

        if not os.path.exists(os.path.join(structure_path, 'structure.pdb')):
            logger.error('Missing protein-ligand complex build.')
            raise FileNotFoundError('Missing protein-ligand complex build.')

        outfile = os.path.join(self.lig_scaff, 'ligand.csv')

        # We only want to exclude sidechain atoms
        selection = 'not segid L and not ('

        for position in self.mutations:
            selection += f"chain {position['chain']} and resid {position['resid']} and not ({backbone_atoms}) or "

        selection = selection[0:-4] + ')'
        struc, prm = load_ff_parameters(structure_path=structure_path, forcefield=self.forcefield)
        ffev = FFEvaluate(struc, prm, betweensets=(selection, 'segid L'))

        self_nrgs = np.zeros((nposes, len(self.energy_terms)))

        with mp.Pool(processes=self.ncpus) as pool:
            with tqdm(total=nposes, desc='Ligand/Scaffold') as pbar:
                for index, energy in enumerate(pool.imap(
                        partial(
                            self.calculate_energy,
                            struc=struc,
                            ffev=ffev,
                            pose_coords=pose_coords
                        ), [(pose,) for pose in range(nposes)],
                        chunksize=calculate_chunks(nposes=nposes, ncpus=self.ncpus))):
                    self_nrgs[index] = energy
                    pbar.update()
        if self.intra:
            # Generate FFEvaluate object for scoring ligand conformation
            struc.filter('segid L', _logger=False)
            ffev = FFEvaluate(struc, prm)

            with mp.Pool(processes=self.ncpus) as pool:
                with tqdm(total=nposes) as pbar:
                    pbar.set_description('Ligand_Internal')
                    for index, energy in enumerate(pool.imap(partial(
                            self.calculate_energy,
                            struc=struc,
                            ffev=ffev,
                            pose_coords=pose_coords), [(pose,) for pose in range(nposes)],
                            chunksize=calculate_chunks(nposes=nposes, ncpus=self.ncpus))):
                        self_nrgs[index] += energy
                        pbar.update()
        # Save data as csv
        write_energies(outpath=outfile,
                       energies=self_nrgs,
                       energy_terms=self.energy_terms,
                       name_a='ligand_pose',
                       nconfs_a=nposes)

    def run_pair_nrg(self, pose_coords: np.ndarray, nposes: int) -> NoReturn:
        """
        Calculate sidechain ligand energies

        Parameters
        ----------
        pose_coords: np.ndarray
            Ligand pose coordinates
        nposes: int
            Number of ligand poses
        """
        from pocketoptimizer.utility.utils import load_ff_parameters, write_energies, calculate_chunks
        from pocketoptimizer.utility.molecule_types import backbone_atoms

        structure_path_prepend = os.path.join(self.work_dir, 'scaffold', self.forcefield, 'protein_params')
        for position in self.mutations:
            chain, resid = position['chain'], position['resid']
            rotamer_path_prepend = os.path.join(self.rotamer_path, f'{chain}_{resid}')
            for resname in position['mutations']:
                output_file = os.path.join(self.lig_side, f'ligand_{chain}_{resid}_{resname}.csv')
                # Check if pairwise-interaction energy file already exists
                if os.path.isfile(output_file):
                    logger.info(f'Ligand-Sidechain interaction energy with residue: {chain}_{resid}_{resname} already computed.')
                    continue
                logger.info(
                    f'Ligand-Sidechain interaction energy with residue: {chain}_{resid}_{resname} not computed yet.')
                # Create Molecule object from rotamers of each residue
                # Check if rotamers are computed for all mutations
                try:
                    rotamer_mol = Molecule(os.path.join(rotamer_path_prepend, f'{resname}.pdb'))
                    rotamer_mol.remove(backbone_atoms, _logger=False)
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

                lig_sidechain = struc.copy()
                lig_sidechain.filter(sel=f'(chain {chain} and resid {resid} and not ({backbone_atoms})) or segid L', _logger=False)

                ffev = FFEvaluate(lig_sidechain, prm,
                                  betweensets=(f'chain {chain} and resid {resid} and not ({backbone_atoms})', 'segid L'))

                batch = list(itertools.product(*[np.arange(nposes), np.arange(nrots)]))
                # Create array (nposes * nrots * 2 (vdw and es))
                pair_nrgs = np.zeros((nposes, nrots * len(self.energy_terms)))

                with tqdm(total=len(batch), desc=f'Ligand/{chain}_{resid}_{resname}') as pbar:
                    with mp.Pool(processes=self.ncpus) as pool:
                        for index, energy in enumerate(pool.imap(
                                partial(
                                    self.calculate_energy,
                                    struc=lig_sidechain,
                                    ffev=ffev,
                                    pose_coords=pose_coords,
                                    rot_coords=rotamer_mol.coords,
                                    chain=chain,
                                    resid=resid
                                ), batch,
                                chunksize=calculate_chunks(nposes=len(batch), ncpus=self.ncpus))):
                            pair_nrgs[batch[index][0], batch[index][1]*2:batch[index][1]*2 + 2] = energy
                            pbar.update()

                # Save data as csv
                write_energies(outpath=output_file,
                               energies=pair_nrgs,
                               energy_terms=self.energy_terms,
                               name_a='ligand_pose',
                               nconfs_a=nposes,
                               name_b=resname,
                               nconfs_b=nrots)

    def run_ff_scorer(self) -> NoReturn:
        """
        Run ligand force field scoring
        """
        logger.info(f'Score ligand interactions using {self.forcefield} force field.')

        os.makedirs(self.lig_scaff, exist_ok=True)
        os.makedirs(self.lig_side, exist_ok=True)

        try:
            ligand_poses = Molecule(self.ligand_poses_pdb)
            ligand_poses.read(self.ligand_poses_xtc)
        except:
            logger.error('Missing ligand poses.')
            raise FileNotFoundError('Missing ligand poses.')

        pose_coords = ligand_poses.coords
        nposes = pose_coords.shape[-1]

        if os.path.isfile(os.path.join(self.lig_scaff, 'ligand.csv')):
            logger.info(f'Ligand-Scaffold/Self interaction energy already computed.')

        else:
            logger.info(f'Ligand-Scaffold/Self interaction energy not computed yet.')
            self.run_self_nrg(pose_coords=pose_coords, nposes=nposes)

        self.run_pair_nrg(pose_coords=pose_coords, nposes=nposes)
        logger.info('Ligand scoring was successful.')
