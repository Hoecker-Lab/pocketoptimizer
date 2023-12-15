import os
import numpy as np
from typing import NoReturn
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import pandas as pd

from pocketoptimizer.utility.utils import create_pairs, Storer
from pocketoptimizer.utility.index_mapper import IndexMapper


class EnergyReader(Storer):
    """
    Class providing functionality to read self & pair energies from
    "raw" energy files and make the energy values accessible
    """

    def __init__(self, index_mapper: IndexMapper, **kwargs):
        """
        Constructor

        Parameters
        ----------
        index_mapper: :class: IndexMapper object
            IndexMapper defining the current design
        scorer: str
            Scoring method used to score ligand/protein interactions
        """
        super().__init__(**kwargs)
        self.index_mapper = index_mapper

        # {'B_41':np.array([13.234, 4234.42, 42.234,...]),
        #  'ligand':np.array([]), 'water3':np_array[-2.0, 0], ...}
        self._self_energies = {}

        # {('A_57', 'B_41'):np.array([13.234, 4234.42, 42.234,...],
        #  ('A_57, 'ligand'):np.array([23.424, ...}
        self._pair_energies = {}

        self._conformer_counts = {}

    def read_energies(self) -> NoReturn:
        """
        Reads the energies from the raw energy files and stores them
        """
        self.read_self_energies()
        self.read_pairwise_energies()

    def read_self_energies(self) -> NoReturn:
        """
        Read the self-interaction energies for ligand conformer/sidechain rotamer with the fixed scaffold part and store
        in numpy arrays as values in self._self_energies dictionary, where the keys are positions ('A_44', 'ligand',..)
        Stack along axis 1 if position already exists in dictionary
        """
        for position in self.mutations:
            chain, resid = position['chain'], position['resid']
            pose = f'{chain}_{resid}'
            # Create 1D array of zeros from number of conformations for position
            nconfs = self.index_mapper.get_conf_count_for_pos(pose)
            self._self_energies[pose] = np.zeros(nconfs)
            i = 0

            for resname in position['mutations']:
                # Filename of csv sdiechain files
                csv = os.path.join(self.side_scaff, f'{chain}_{resid}_{resname}.csv')
                # Read in csv as pandas dataframe
                df = pd.read_csv(csv, delimiter='\t', index_col=0)
                # create numpy array from dataframe summing up energies
                # get energy offset
                nrg_offset = len(df.columns)
                energies = np.transpose(np.array(df.groupby((np.arange(len(df.columns)) // nrg_offset) + 1, axis=1).sum()))
                # Fill array with energies
                self._self_energies[pose][i:i+energies.shape[-1]] = energies[0]
                i += energies.shape[-1]

        nconfs = self.index_mapper.get_conf_count_for_pos('ligand')
        self._self_energies['ligand'] = np.zeros(nconfs)
        # Filename of csv ligand file
        csv = os.path.join(self.lig_scaff, 'ligand.csv')
        # Read in csv as pandas dataframe
        df = pd.read_csv(csv, delimiter='\t', index_col=0)
        # get energy offset depending on scoring function
        nrg_offset = len(df.columns)
        # create numpy array from dataframe summing up energies of respective scoring function
        energies = np.transpose(np.array(df.groupby((np.arange(len(df.columns)) // nrg_offset) + 1, axis=1).sum()))
        self._self_energies['ligand'][0:energies.shape[-1]] = energies[0]

    def read_pairwise_energies(self) -> NoReturn:
        """
        Read pairwise-interaction energies for sidechain rotamer/sidechain rotamer (Only vdw and electrostatic energies).
        Read pairwise-interaction energies for sidechain rotamer/ligand conformer (N energy terms as calculated by the respective scoring function)
        The energies are written per position and stored in an numpy array which has the dimensions
        NxM, where N is the number of rotamers for all residues at position A and M is the number of rotamers/ligand
        conformers for all residues at position B
        """
        sorted_pair_permutations = create_pairs(self.mutations)

        # Read side chain/side chain energies
        for mutation_a, mutation_b in sorted_pair_permutations:
            # Get name of each residue from each pair
            chain_a, resid_a = mutation_a.split('_')[0:2]
            chain_b, resid_b = mutation_b.split('_')[0:2]
            resname_a = mutation_a.split('_')[2]
            resname_b = mutation_b.split('_')[2]
            pos_a = f'{chain_a}_{resid_a}'
            pos_b = f'{chain_b}_{resid_b}'

            # Check if array already exists for pose
            if not (pos_a, pos_b) in self._pair_energies:
                # Create 2D array of zeros in 1st dimension number of conformations for pose_a and in 2nd dimension number of conformations for pose_b
                nconfs_a, nconfs_b = self.index_mapper.get_conf_count_for_pos(pos_a), self.index_mapper.get_conf_count_for_pos(pos_b)
                # Create array for energy values of the two positions
                self._pair_energies[(pos_a, pos_b)] = np.zeros((nconfs_a, nconfs_b))
                # start new row and column counters
                row = 0
                column = 0
                tmp_resname_a = resname_a
            csv = os.path.join(self.side_side,
                               f'{chain_a}_{resid_a}_{resname_a}_{chain_b}_{resid_b}_{resname_b}.csv')
            # Read in csv as pandas dataframe
            df = pd.read_csv(csv, delimiter='\t', index_col=0)
            # create numpy array from dataframe summing up vdw and es energies
            energies = np.array(df.groupby((np.arange(len(df.columns)) // 2) + 1, axis=1).sum())
            # resname_b has changed, go to next columns of array and increase column counter, save row to jump to next row, if resname_a changes
            if resname_a == tmp_resname_a:
                self._pair_energies[(pos_a, pos_b)][row:row+self.index_mapper.get_conf_count_for_res(pos_a, resname_a), column:column+energies.shape[-1]] = energies
                column += energies.shape[-1]
                old_row = self.index_mapper.get_conf_count_for_res(pos_a, resname_a)
            # If resname_a has changed, go to next row and initialize new column counter at 0
            else:
                row += old_row
                self._pair_energies[(pos_a, pos_b)][row:row+self.index_mapper.get_conf_count_for_res(pos_a, resname_a), 0:energies.shape[-1]] = energies
                column = energies.shape[-1]
                old_row = self.index_mapper.get_conf_count_for_res(pos_a, resname_a)
                tmp_resname_a = resname_a

        # Read ligand/side chain energies
        for position in self.mutations:
            pos_a = f'{position["chain"]}_{position["resid"]}'
            pos_b = 'ligand'
            # Create 2D array of zeros in 1st dimension number of conformations for pose_a and in 2nd dimension number of conformations for pose_b
            nconfs_a, nconfs_b = self.index_mapper.get_conf_count_for_pos(pos_a), self.index_mapper.get_conf_count_for_pos(pos_b)
            # Create array for energy values of the two positions
            self._pair_energies[(pos_a, pos_b)] = np.zeros((nconfs_a, nconfs_b))
            row = 0
            for resname in position['mutations']:
                # Filename of csv ligand file
                csv = os.path.join(self.lig_side, f'ligand_{position["chain"]}_{position["resid"]}_{resname}.csv')
                # Read in csv as pandas dataframe
                df = pd.read_csv(csv, delimiter='\t', index_col=0)
                nrg_offset = int(len(df.columns)/(int(df.columns[-1].split('_')[1]) + 1))
                # create numpy array from dataframe summing up vdw and es energies
                energies = np.transpose(np.array(df.groupby((np.arange(len(df.columns)) // nrg_offset) + 1, axis=1).sum()))
                # Fill array with energies
                self._pair_energies[(pos_a, pos_b)][row:row+self.index_mapper.get_conf_count_for_res(pos_a, resname), 0:energies.shape[-1]] = energies
                row += self.index_mapper.get_conf_count_for_res(pos_a, resname)

    def get_self_energy(self, pos: str, index: int) -> float:
        """
        Get a self energy value for a position and residue

        Parameters
        ----------
        pos: str
            Design position name
        index: int
            Conformer index of residue

        Returns
        -------
        Self-interaction energy value for conformation at position
        """
        return self._self_energies[pos][index] * self.index_mapper.get_scaling_factor(pos)

    def get_pair_energy(self, pos1: str, pos2: str, i1: int, i2: int) -> float:
        """
        Get a pair energy value

        Parameters
        ----------
        pos1: str
            Name of design position 1
        pos2: str
            Name of design position 2
        i1: int
            Conformer index of res at position 1
        i2: int
            Conformer index of res at position 2

        Returns
        -------
        Pairwise-interaction energy value for two conformations at two positions
        """
        sf1 = self.index_mapper.get_scaling_factor(pos1)
        sf2 = self.index_mapper.get_scaling_factor(pos2)
        return self._pair_energies[(pos1, pos2)][i1][i2] * sf1 * sf2

    def get_max_self_energy(self, pos: str) -> float:
        """
        Get the maximal self energy at a position

        Parameters
        ----------
        pos: str
            Design position name

        Returns
        -------
        Maximum self-interaction energy value for position
        """
        return max(self._self_energies[pos])
