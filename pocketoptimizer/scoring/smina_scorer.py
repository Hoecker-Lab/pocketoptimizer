import os
import shutil
import subprocess
import tempfile as tf
from typing import List, Dict, NoReturn
import logging
import numpy as np
from moleculekit.molecule import Molecule
from tqdm.auto import tqdm

from pocketoptimizer.utility.utils import Storer

logger = logging.getLogger(__name__)


class SminaScorer(Storer):
    """
    SminaScorer class for protein/ligand interactions
    Read scaffold
    Remove side chains at design positions from scaffold
    Write pdb
    Combine ligand mol2 files into one file and multiprocess scoring with smina
    Parse output
    Write output into .csv files

    """
    # Actually there should be conversion of receptor to pdbqt format
    # But since it does not make any difference in scores it has been removed

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

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
        if self.scorer not in self.smina_weights.keys():
            logger.error(f'{self.scorer} is not a valid scoring method. Valid methods are:'
                         f' vina, vinardo, ad4_scoring.')
            raise ValueError(f'{self.scorer} is not a valid scoring method. Valid methods are:'
                             f' vina, vinardo, ad4_scoring.')
        self.terms = [weight for weight in self.smina_weights[self.scorer].keys()]
        self.weights = [weight for weight in self.smina_weights[self.scorer].values()]

    def parse_smina(self, smina_output: str, nposes: int) -> np.ndarray:
        """
        Smina output parsing function

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
        nrgs = np.zeros(len(self.weights))
        scores = np.zeros((nposes, len(self.smina_weights[self.scorer])))

        j = 0
        for line in smina_output.split('\n'):
            # Read out the energies for each pose
            if line.startswith('## ') and not line.startswith('## Name'):
                for i, nrg in enumerate(line.strip().split()[2:]):
                    nrgs[i] = float(nrg) * self.weights[i]
                scores[j] = nrgs
                j += 1
        return scores

    def prepare_scaffold(self) -> NoReturn:
        """
        Prepares a scaffold where all side chains of mutation positions are removed

        Parameters
        ----------
        mutations: List[Dict]
            List of dictionaries defining the mutations at different positions
        """
        from pocketoptimizer.utility.molecule_types import backbone_atoms

        scaffold = Molecule(self.built_scaffold)
        for position in self.mutations:
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

    def score_smina(self, receptor: str, ligand: str, nposes: int) -> np.ndarray:
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

        Returns
        -------
        Array containing weighted energy components for all poses
        """

        score_command = [
            self.smina,
            '--receptor', receptor,
            '--ligand', ligand,
            '--score_only',
            '--cpu', str(self.ncpus),
            '--scoring', self.scorer
        ]
        process = subprocess.run(score_command, capture_output=True)
        return self.parse_smina(process.stdout.decode('ascii'), nposes=nposes)

    def run_smina_scorer(self, _keep_tmp: bool = False) -> NoReturn:
        """
        Wrapper to run smina ligand protein interaction energy scoring
        Reading in ligand structure (ligand.mol2) and ligand poses from starting pose ligand_poses.pdb and ligand_poses.xtc trajectory

        Parameters
        ----------
        _keep_tmp: bool
            If False deleting the temporary directory afterwards [default: False]
        """
        from pocketoptimizer.utility.utils import write_energies
        logger.info(f'Score ligand interactions using {self.scorer}.')

        os.makedirs(self.lig_scaff, exist_ok=True)
        os.makedirs(self.lig_side, exist_ok=True)

        self.tmp_dir = tf.mkdtemp(dir=self.tmp_dir, prefix=f'calculateLigand{self.scorer}_')
        os.chdir(self.tmp_dir)

        ligand = Molecule(self.ligand_protonated)
        try:
            ligand_poses = Molecule(self.ligand_poses_pdb)
            ligand_poses.read(self.ligand_poses_xtc)
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
            self.prepare_scaffold()
            with tqdm(total=1, desc='Ligand/Scaffold') as pbar:
                self_nrgs = self.score_smina(receptor='pocketless.pdb',
                                             ligand=lig_outfile,
                                             nposes=nposes)
                pbar.update()

            # Save data as csv
            energy_terms = self.terms

            write_energies(outpath=lig_scaffold_outfile,
                               energies=self_nrgs,
                               energy_terms=energy_terms,
                               name_a='ligand_pose',
                               nconfs_a=nposes)

        # Set the torsion weight factor to 0, in order not to count it for every pairwise interaction
        self.smina_weights['ad4_scoring']['num_tors_add'] = 0.0

        # Score ligand against sidechains
        # Loop over all mutations/flexible residues
        for position in self.mutations:
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
                                                    nposes=nposes)
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
