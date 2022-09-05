import os
import shutil
import subprocess
import tempfile as tf
from typing import List, Dict, Union, NoReturn
import logging
import numpy as np
from moleculekit.molecule import Molecule
from tqdm.auto import tqdm

logger = logging.getLogger(__name__)


class SminaScorer:
    """
    SminaScorer class for ligand protein interactions (ligand sidechain and ligand scaffold)
    Read scaffold
    Remove side chains at design positions from scaffold
    Write pdb
    Combine ligand mol2 files into one file and multiprocess scoring with smina
    Parse output
    Write output into .csv files

    """
    # Actually there should be conversion of receptor to pdbqt format
    # But since it does not make any difference in scores it has been removed

    def __init__(self, ligand: str, ligand_poses: str, scaffold: str, project_path: str, scoring_method: str,
                 forcefield: str, rotamer_path: str, smina_path: str, tmp_dir: str = '/tmp/'):
        """
        Initializes a :class: SminaScorer object

        Parameters
        ----------
        ligand: str
            Path to ligand.mol2 file with sybyl atom types
        ligand_poses: str
            Path to ligand_poses.pdb file
        scaffold: str
            Path to protein scaffold file
        pocketless_scaffold: str
            Path to Structure without pocket sidechains
        project_path: str
            Path of the whole project
        scoring_method: str
            Scoring method to use
        forcefield: str
            Force field used to parameterize the ligand
        rotamer_path: str
            Path to sampled side chain rotamers
        smina_path: str
            Path to smina binary
        tmp_dir: str
            path to temporary directory [default: '/tmp/']
        """
        self.ligand = ligand
        self.ligand_poses = ligand_poses
        self.scaffold = scaffold
        self.pocketless_scaffold = 'pocketless.pdb'
        self.project_path = project_path
        self.lig_scaff = ''
        self.lig_side = ''
        self.scoring_method = scoring_method

        # The weights of the smina energy components as provided by the software
        self.smina_weights = {
            'vina': {
                'gauss1': -0.035579,
                'gauss2': -0.005156,
                'repulsion': 0.840245,
                'hydrophobic': -0.035069,
                'non_dir_h_bond': -0.587439,
                'num_tors_div': 1.923
            },
            'vinardo': {
                'gauss': -0.045,
                'repulsion': 0.8,
                'hydrophobic': -0.035,
                'non_dir_h_bond': -0.6,
                'num_tors_div': 0.0
            },
            'ad4_scoring': {
                'vdw': 0.156,
                'non_dir_h_bond_lj': -0.0974,
                'ad4_solvation': 0.1159,
                'electrostatic': 0.1465,
                'num_tors_add': 0.2744
            }
        }
        self.terms = [weight for weight in self.smina_weights[self.scoring_method].keys()]
        self.weights = [weight for weight in self.smina_weights[self.scoring_method].values()]
        if self.scoring_method not in self.smina_weights.keys():
            logger.error(f'{self.scoring_method} is not a valid scoring method. Valid methods are:'
                         f' vina, vinardo, ad4_scoring.')
            raise ValueError(f'{self.scoring_method} is not a valid scoring method. Valid methods are:'
                             f' vina, vinardo, ad4_scoring.')
        self.forcefield = forcefield
        self.rotamer_path = rotamer_path
        self.smina = smina_path
        self.tmp_dir = tmp_dir

    def parse_smina(self, smina_output: str, nposes: int, intra: bool = False) -> np.ndarray:
        """
        Smina output parsing function

        Parameters
        ----------
        smina_output: str
            Standard output captured from calling smina
        nposes: int
            Number of ligand poses
        intra: bool
            Whether to parse internamolecular energy

        Returns
        -------
        Array containing weighted energy components for all poses
        """
        if intra:
            nrgs = np.zeros(len(self.weights) + 1)
            scores = np.zeros((nposes, len(self.smina_weights[self.scoring_method]) + 1))
        else:
            nrgs = np.zeros(len(self.weights))
            scores = np.zeros((nposes, len(self.smina_weights[self.scoring_method])))

        j = 0
        for line in smina_output.split('\n'):
            # Read out the energies for each pose
            if intra:
                if line.startswith('Intramolecular energy:'):
                    nrgs[0] = float(line.strip().split()[2])
                elif line.startswith('## ') and not line.startswith('## Name'):
                    for i, nrg in enumerate(line.strip().split()[2:]):
                        nrgs[i + 1] = float(nrg) * self.weights[i]
                    scores[j] = nrgs
                    j += 1
            else:
                if line.startswith('## ') and not line.startswith('## Name'):
                    for i, nrg in enumerate(line.strip().split()[2:]):
                        nrgs[i] = float(nrg) * self.weights[i]
                    scores[j] = nrgs
                    j += 1
        return scores

    def prepare_scaffold(self, mutations: List[Dict[str, Union[str, List[str]]]]) -> NoReturn:
        """
        Prepares a scaffold where all side chains of mutation positions are removed

        Parameters
        ----------
        mutations: List[Dict]
            List of dictionaries defining the mutations at different positions
        """
        from pocketoptimizer.utility.molecule_types import backbone_atoms

        scaffold = Molecule(self.scaffold)
        for position in mutations:
            scaffold.remove(f"chain {position['chain']} and resid {position['resid']} and not ({backbone_atoms})", _logger=False)
        scaffold.write('tmp.pdb')
        # Remove END from pdb file because babel sucks
        with open('tmp.pdb', 'r') as infile, open('pocketless.pdb', 'w') as outfile:
            lines = infile.readlines()
            for index, line in enumerate(lines):
                if line.startswith('CONECT'):
                    break
            outfile.writelines(lines[:index])
        assert os.path.isfile('pocketless.pdb')

    def prepare_sidechain(self, rotamers: Molecule, chain: str, resid: str, resname: str) -> NoReturn:
        """
        Write pdb files for each rotamer

        Parameters
        ----------
        rotamers: :class: moleculekit.Molecule object
            Containing all rotamers for a residue
        chain: str
            Chain identifier of residue
        resid: str
            Residue identifier of residue
        resname: str
            Name of residue
        """
        for rotamer in range(rotamers.coords.shape[-1]):
            res_pose_name = f'{chain}_{resid}_{resname}_{rotamer}.pdb'
            if not os.path.isfile(res_pose_name):
                rot = rotamers.copy()
                # Needs shape (natoms, 3, 1)
                rot.coords = rotamers.coords[:, :, rotamer].reshape(rotamers.coords.shape[0], 3, 1)
                rot.write(res_pose_name)

    def combine_ligand(self, poses: int, outfile: str) -> NoReturn:
        """
        Combines ligand_pose.mol2 files into single file to minimize read/write accesses when scoring with smina/gnina/oddt

        poses: int
            Number of poses to combine
        outfile: str
            Name of combined output file
        """
        # Just parse the mol2 files together
        with open(outfile, 'w') as combined_file:
            for pose_id in range(poses):
                with open(f'ligand_pose_{pose_id}.mol2', 'r') as pose_file:
                    for line in pose_file:
                        combined_file.write(line)

    def score_smina(self, receptor: str, ligand: str, nposes: int, intra: bool = False, ncpus: int = 1) -> Dict[str, float]:
        """
        Runs Smina scoring process for a receptor and ligand/s

        Parameters
        ----------
        receptor: str
            Path to input protein structure (either side chain removed scaffold or sidechain rotamer)
        ligand: str
            Ligand structure/s
        nposes: int
            Number of ligand poses
        intra: bool
            Whether to parse internamolecular energy
        ncpus: int
            Number of CPUs to use for scoring [default: 1]

        Returns
        -------
        List of dictionaries containing different energy terms as keys and weighted energy values as values ({'gauss1': energy,..})
        """

        score_command = [
            self.smina,
            '--receptor', receptor,
            '--ligand', ligand,
            '--score_only',
            '--cpu', str(ncpus),
            '--scoring', self.scoring_method
        ]
        process = subprocess.run(score_command, capture_output=True)
        return self.parse_smina(process.stdout.decode('ascii'), nposes=nposes, intra=intra)

    def run_smina_scorer(self, mutations: List[Dict[str, Union[str, List[str]]]],
                         ncpus: int = 1, _keep_tmp: bool = False) -> None or NoReturn:
        """
        Wrapper to run smina ligand protein interaction energy scoring
        Reading in ligand structure (ligand.mol2) and ligand poses from starting pose ligand_poses.pdb and ligand_poses.xtc trajectory

        Parameters
        ----------
        mutations: List[Dict]
            List of dictionaries defining the mutations at different positions
        ncpus: int
            Number of CPUs to use for scoring [default: number of CPUs in the system]
        _keep_tmp: bool
            If False deleting the temporary directory afterwards [default: False]
        """
        from pocketoptimizer.utility.utils import write_energies

        logger.info(f'Score ligand interactions using {self.scoring_method}.')

        # Output ligand scaffold energy
        self.lig_scaff = os.path.join(self.project_path,
                                      'energies',
                                      f'{self.forcefield}_{self.rotamer_path.split("/")[-1]}',
                                      f'ligand_scaffold_{self.scoring_method}')
        os.makedirs(self.lig_scaff, exist_ok=True)
        # Output ligand side chain energies
        self.lig_side = os.path.join(self.project_path,
                                     'energies',
                                     f'{self.forcefield}_{self.rotamer_path.split("/")[-1]}',
                                     f'ligand_sidechain_{self.scoring_method}')
        os.makedirs(self.lig_side, exist_ok=True)

        self.tmp_dir = tf.mkdtemp(dir=self.tmp_dir, prefix=f'calculateLigand{self.scoring_method}_')
        os.chdir(self.tmp_dir)

        ligand = Molecule(self.ligand)

        try:
            ligand_poses = Molecule(self.ligand_poses)
            ligand_poses.read(f'{os.path.splitext(self.ligand_poses)[0]}.xtc')
        except:
            logger.error('Missing ligand poses.')
            raise FileNotFoundError('Missing ligand poses.')

        nposes = ligand_poses.coords.shape[-1]

        # Create combined ligand pose file
        for pose in range(nposes):
            ligand.coords = ligand_poses.coords[:, :, pose].reshape(ligand_poses.coords.shape[0], 3, 1)
            pose = f'ligand_pose_{pose}.mol2'
            ligand.write(pose)
        # Combine poses
        lig_outfile = os.path.join(self.tmp_dir, f'ligand_poses.mol2')
        self.combine_ligand(poses=nposes,
                            outfile=lig_outfile)

        # Score ligand against fixed scaffold
        lig_scaffold_outfile = os.path.join(self.lig_scaff, 'ligand.csv')
        if os.path.isfile(lig_scaffold_outfile):
            logger.info(f'Ligand-Scaffold/Self interaction energy already computed.')
        else:
            logger.info(f'Ligand-Scaffold/Self interaction energy not computed yet.')
            logger.info('Prepare fixed scaffold.')
            self.prepare_scaffold(mutations=mutations)
            with tqdm(total=1, desc='Ligand/Scaffold') as pbar:
                self_nrgs = self.score_smina(receptor=self.pocketless_scaffold,
                                             ligand=lig_outfile,
                                             nposes=nposes,
                                             intra=False,
                                             ncpus=ncpus)
                pbar.update()
            # Save data as csv
            write_energies(outpath=lig_scaffold_outfile,
                           energies=self_nrgs,
                           #energy_terms=['Intra'] + self.terms,
                           energy_terms=self.terms,
                           name_a='ligand_pose',
                           nconfs_a=nposes)

        # Set the torsion weight factor to 0, in order not to count it for every pairwise interaction
        self.smina_weights['ad4_scoring']['num_tors_add'] = 0.0

        # Score ligand against sidechains
        # Loop over all mutations/flexible residues
        for position in mutations:
            chain, resid = position['chain'], position['resid']
            for resname in position['mutations']:
                lig_sidechain_outfile = os.path.join(f"{self.lig_side}", f"ligand_{chain}_{resid}_{resname}.csv")
                # Check if pairwise-interaction energy file already exists
                if os.path.isfile(lig_sidechain_outfile):
                    logger.info(f'Ligand-Sidechain interaction energy for residue: {chain}_{resid}_{resname} already computed.')
                    continue
                logger.info(
                    f'Ligand-Sidechain interaction energy for residue: {chain}_{resid}_{resname} not computed yet.')
                # Create Molecule object from rotamers of each residue
                # Check if rotamers are computed for all mutations
                try:
                    sidechain = Molecule(os.path.join(self.rotamer_path, f"{chain}_{resid}", f"{resname}.pdb"))

                except FileNotFoundError:
                    logger.error(f'Missing rotamer for residue: {chain}_{resid}_{resname}.')
                    raise FileNotFoundError(f'Missing rotamer for residue: {chain}_{resid}_{resname}.')

                # select only sidechain
                sidechain.filter('sidechain', _logger=False)

                self.prepare_sidechain(sidechain, chain, resid, resname)

                # Number of rotamers to score
                nrots = sidechain.coords.shape[-1]
                # Number of terms in the scoring function
                nterms = len(self.terms)

                # Create array (nposes * nrots * Number of energy terms)
                pair_nrgs = np.zeros((nposes, nrots * nterms))

                logger.info(f'Loop over rotamers of residue: {chain}_{resid}_{resname}.')
                with tqdm(total=nrots, desc=f'Ligand/{chain}_{resid}_{resname}') as pbar:
                    for rot in np.arange(nrots):
                        res_pdb = f'{chain}_{resid}_{resname}_{rot}.pdb'
                        pair_nrg = self.score_smina(receptor=res_pdb,
                                                    ligand=lig_outfile,
                                                    nposes=nposes,
                                                    intra=False,
                                                    ncpus=ncpus)
                        pair_nrgs[:, rot*nterms:rot*nterms + nterms] = pair_nrg
                        pbar.update()

                # Save data as csv
                write_energies(outpath=lig_sidechain_outfile,
                               energies=pair_nrgs,
                               energy_terms=self.terms,
                               name_a='ligand_pose',
                               nconfs_a=nposes,
                               name_b=resname,
                               nconfs_b=nrots)

        if not _keep_tmp:
            if os.path.isdir(self.tmp_dir):
                shutil.rmtree(self.tmp_dir)

        logger.info('Ligand scoring was successful.')
