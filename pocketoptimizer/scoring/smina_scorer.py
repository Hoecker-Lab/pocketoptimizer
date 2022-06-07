import multiprocessing as mp
import os
import shutil
import subprocess
import tempfile as tf
from typing import List, Dict, Union, NoReturn
import logging
import numpy as np
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

logger = logging.getLogger('pocketoptimizer.scoring.smina_scorer')


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
                 forcefield: str, rotamer_path: str, smina_path: str, obabel_path: str):
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
        obabel_path: str
            Path to obabel binary
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
            'dkoes_scoring': {
                'vdw': 0.0099,
                'non_dir_h_bond': -0.153055,
                'ad4_solvation': 0.048934,
                'num_tors_sqr': 0.317267,
                'constant_term': -2.46902
            },
            'ad4_scoring': {
                'vdw': 0.156,
                'non_dir_h_bond_lj': -0.0974,
                'ad4_solvation': 0.1159,
                'electrostatic': 0.1465,
                'num_tors_add': 0.2744
            }
        }
        self.forcefield = forcefield
        self.rotamer_path = rotamer_path
        self.smina = smina_path
        self.obabel = obabel_path
        if self.scoring_method not in self.smina_weights.keys():
            logger.error(f'{self.scoring_method} is not a valid scoring method. Valid methods are:'
                             f' vina, vinardo, dkoes_scoring, ad4_scoring.')
            raise ValueError(f'{self.scoring_method} is not a valid scoring method. Valid methods are:'
                             f' vina, vinardo, dkoes_scoring, ad4_scoring.')

    def parse_smina(self, smina_output: str, nposes: int) -> np.ndarray:
        """
        Parses smina affinity values

        #################################################################\n
         # If you used AutoDock Vina in your work, please cite:          #\n
         #                                                               #\n
         # O. Trott, A. J. Olson,                                        #\n
         # AutoDock Vina: improving the speed and accuracy of docking    #\n
         # with a new scoring function, efficient optimization and       #\n
         # multithreading, Journal of Computational Chemistry 31 (2010)  #\n
         # 455-461                                                       #\n
         #                                                               #\n
         # DOI 10.1002/jcc.21334                                         #\n
         #                                                               #\n
         # Please see http://vina.scripps.edu for more information.      #\n
         #################################################################\n
         Detected 8 CPUs\n
         Reading input ... done.\n
         Setting up the scoring function ... done.\n
         Affinity: -0.60387 (kcal/mol)\n
         Intermolecular contributions to the terms, before weighting:\n
             gauss 1     : 65.69073\n
             gauss 2     : 633.94241\n
             repulsion   : 8.09557\n
             hydrophobic : 6.89632\n
             Hydrogen    : 2.71140\n
        Parameters
        ----------
        smina_output: str
            Standard output captured from calling smina
        nposes: int
            Number of ligand poses

        Returns
        -------
        Array containing weighted energy components for all poses
        """
        weights = [weight for weight in self.smina_weights[self.scoring_method].values()]
        scores = np.zeros((nposes, len(self.smina_weights[self.scoring_method])))

        j = 0
        for line in smina_output.split('\n'):
            # Read out the energies for each pose
            if line.startswith('## ') and not line.startswith('## Name'):
                nrgs = np.zeros(len(weights))
                for i, nrg in enumerate(line.strip().split()[2:]):
                    nrgs[i] = float(nrg) * weights[i]
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
        scaffold = Molecule(self.scaffold)
        for position in mutations:
            scaffold.remove(f"chain {position['chain']} and resid {position['resid']} and sidechain", _logger=False)
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
        with open(outfile, 'w') as combined_file:
            for i in range(poses):
                with open(f'ligand_pose_{i}.mol2', 'r') as pose_file:
                    for line in pose_file:
                        combined_file.write(line)

    def score_smina(self, receptor: str, ligand: str, nposes: int, ncpus: int = 1) -> Dict[str, float]:
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
        ncpus: int
            Number of CPUs to use for scoring [default: 1]

        Returns
        -------
        List of dictionaries containing different energy terms as keys and weighted energy values as values ({'gauss1': energy,..})
        """

        command = [
            self.smina,
            '--receptor', receptor,
            '--ligand', ligand,
            '--score_only',
            '--cpu', str(ncpus),
            '--scoring', self.scoring_method
        ]
        process = subprocess.run(command, capture_output=True)
        return self.parse_smina(process.stdout.decode('ascii'), nposes=nposes)

    def run_smina_scorer(self, mutations: List[Dict[str, Union[str, List[str]]]], temp_dir: str = '/tmp/', ncpus: int = mp.cpu_count(), _keep_tmp: bool = False) -> None or NoReturn:
        """
        Wrapper to run smina ligand protein interaction energy scoring
        Reading in ligand structure (ligand.mol2) and ligand poses from starting pose ligand_poses.pdb and ligand_poses.xtc trajectory

        Parameters
        ----------
        mutations: List[Dict]
            List of dictionaries defining the mutations at different positions
        temp_dir: str
            Name of temporary directory to create a temporary directory in containing
            all .pdb files for rotamers and ligand pose .mol2 files [default: '/tmp/']
        ncpus: int
            Number of CPUs to use for scoring [default: number of CPUs in the system]
        _keep_tmp: bool
            If False deleting the temporary directory afterwards [default: False]
        """
        from pocketoptimizer.utility.utils import write_energies

        logger.info(f'Score ligand interactions using {self.scoring_method}.')
        weights = [weight for weight in self.smina_weights[self.scoring_method].keys()]

        # Output ligand scaffold energy
        self.lig_scaff = os.path.join(self.project_path, 'energies', f'{self.forcefield}_{self.rotamer_path.split("/")[-1]}', f'ligand_scaffold_{self.scoring_method}')
        os.makedirs(self.lig_scaff, exist_ok=True)
        # Output ligand side chain energies
        self.lig_side = os.path.join(self.project_path, 'energies', f'{self.forcefield}_{self.rotamer_path.split("/")[-1]}', f'ligand_sidechain_{self.scoring_method}')
        os.makedirs(self.lig_side, exist_ok=True)

        tmp_dir = tf.mkdtemp(dir=temp_dir, prefix=f'calculateLigand{self.scoring_method}_')
        os.chdir(tmp_dir)

        ligand = Molecule(self.ligand)

        try:
            ligand_poses = Molecule(self.ligand_poses)
            ligand_poses.read(f'{os.path.splitext(self.ligand_poses)[0]}.xtc')
        except:
            logger.error('Missing ligand poses.')
            raise FileNotFoundError('Missing ligand poses.')

        nposes = ligand_poses.coords.shape[-1]

        # Create ligand conformer .mol2 files
        for pose in range(nposes):
            ligand.coords = ligand_poses.coords[:, :, pose].reshape(ligand_poses.coords.shape[0], 3, 1)
            pose_mol2 = f'ligand_pose_{pose}.mol2'
            ligand.write(pose_mol2)
        # Combine mol2 files
        lig_outfile = 'ligand_poses.mol2'
        self.combine_ligand(nposes, lig_outfile)

        # Score ligand against fixed scaffold
        lig_scaffold_outfile = os.path.join(self.lig_scaff, 'ligand.csv')
        if os.path.isfile(lig_scaffold_outfile):
            logger.info(f'Ligand-Scaffold-Interaction-Energies already computed.')
        else:
            logger.info(f'Ligand-Scaffold-Interaction-Energies not computed yet.')
            logger.info('Prepare fixed scaffold.')
            self.prepare_scaffold(mutations=mutations)
            with tqdm(total=1, desc='Ligand/Scaffold') as pbar:
                self_nrgs = self.score_smina(receptor=self.pocketless_scaffold, ligand=lig_outfile, nposes=nposes, ncpus=ncpus)
                pbar.update()

            # Save data as csv
            write_energies(outpath=lig_scaffold_outfile,
                           energies=self_nrgs,
                           energy_terms=weights,
                           name_a='ligand_pose',
                           nconfs_a=nposes
                           )

        # Score ligand against sidechains
        # Loop over all mutations/flexible residues
        for position in mutations:
            chain, resid = position['chain'], position['resid']
            for resname in position['mutations']:
                lig_sidechain_outfile = os.path.join(f"{self.lig_side}", f"ligand_{chain}_{resid}_{resname}.csv")
                # Check if pairwise-interaction energy file already exists
                if os.path.isfile(lig_sidechain_outfile):
                    logger.info(f'Ligand-Sidechain-Interaction-Energies with Residue: {chain}_{resid}_{resname} already computed.')
                    continue
                logger.info(
                    f'Ligand-Sidechain-Interaction-Energies with Residue: {chain}_{resid}_{resname} not computed yet.')
                # Create Molecule object from rotamers of each residue
                # Check if rotamers are computed for all mutations
                try:
                    sidechain = Molecule(os.path.join(self.rotamer_path, f"{chain}_{resid}", f"{resname}.pdb"))

                except FileNotFoundError:
                    logger.error(f'Missing rotamer for residue: {chain}_{resid}_{resname}.')
                    raise FileNotFoundError(
                        f'Missing rotamer for residue: {chain}_{resid}_{resname}.')

                # select only sidechain
                sidechain.filter('sidechain', _logger=False)

                self.prepare_sidechain(sidechain, chain, resid, resname)

                # Number of rotamers to score
                nrots = sidechain.coords.shape[-1]
                # Number of terms in the scoring function
                nterms = len(weights)

                # Create array (nposes * nrots * Number of energy terms)
                pair_nrgs = np.zeros((nposes, nrots * nterms))

                logger.info(f'Loop over Rotamers of Residue: {chain}_{resid}_{resname}.')
                with tqdm(total=nrots, desc=f'Ligand/{chain}_{resid}_{resname}') as pbar:
                    for rot in np.arange(nrots):
                        res_pdb = f'{chain}_{resid}_{resname}_{rot}.pdb'
                        pair_nrg = self.score_smina(receptor=res_pdb, ligand=lig_outfile, nposes=nposes, ncpus=ncpus)
                        pair_nrgs[:, rot*nterms:rot*nterms + nterms] = pair_nrg
                        pbar.update()

                # Save data as csv
                write_energies(outpath=lig_sidechain_outfile,
                               energies=pair_nrgs,
                               energy_terms=weights,
                               name_a='ligand_pose',
                               nconfs_a=nposes,
                               name_b=resname,
                               nconfs_b=nrots)

        if not _keep_tmp:
            if os.path.isdir(tmp_dir):
                shutil.rmtree(tmp_dir)

        logger.info('Ligand scoring was successful.')
