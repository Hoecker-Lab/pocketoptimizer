import os
from typing import List, Dict, Union, NoReturn
import numpy as np
import pandas as pd
from moleculekit.molecule import Molecule

from pocketoptimizer.utility.index_mapper import IndexMapper


class DesignSolution:
    """
    Class that can be used to parse a Sontag solver outfile, extract
       the solutions from it, and generate readable reports and
       structure files of the proposed designs.
    """

    def __init__(self, work_dir: str, mutations: List[Dict[str, Union[str, List[str]]]],
                 forcefield: str, scaffold: str, ligand: str, ligand_poses: str, rotamer_path: str, scorer: str,
                 solution_file: str, index_file: str, energy_file: str):
        """
        Constructor

        Parameters
        ----------
        work_dir: str
            Path to working directory
        mutations: list
            List of all mutations, contains dictionaries as list entries
        forcefield: str
            Forcefield used for energy computations
        scaffold: str
            Path to minimized scaffold structure
        ligand: str
            Path to protonated ligand
        ligand_poses: str
            Path to ligand poses
        rotamer_path: str
            Path to sampeld side chain rotamers
        scorer: str
            Scoring method used to score ligand/protein interactions
        solution_file: str
            The path of the output file of the
            kingsford solver containing the solutions
        index_file: str
            Path of the index file mapping kingsford index
            notation to "readable" position and residue names
        energy_file: str
            Path of the solver input energy file
        """

        self.work_dir = work_dir
        self.mutations = mutations
        self.forcefield = forcefield
        self.scaffold = scaffold
        self.ligand = ligand
        self.ligand_poses = ligand_poses
        self.rotamer_path = rotamer_path
        self.scorer = scorer

        self._raw_conf_indices = []
        # a list of indices for each solution extracted directly
        # from the output file. each is the index of the best
        # conformer for the position [[3, 43, 45...],[..]]

        self._solution_confs = []
        # pattern: [{'B_41':('ARG', 23),...}, ...]
        # for each solution one dict

        self._pairwise_energies = []
        # pattern: [{(3, 5): 10.23, ...},...] . pos indices start at zero
        # for each solution one dict

        self._self_energies = []
        # pattern: [[-123.4, -1.2, 0.0, ...], [...], ...]
        # for each solution one list, indices from index file

        self._detailed_self_energies = []
        # [{'B_41':[-107.003, 0.73742, ...], ...}, ...]
        # energy values in the order as in self._self_energy_components

        self._detailed_pair_energies = []
        # [{('A_32','B_41'):[23.3, ...], \
        #   ('B_41', 'ligand'):[-107.003, 0.73742, ...], ...}, ...]
        # energy values in the order as in self._self_energy_components

        self._self_energy_components = {}
        # maps the self energy of a position to the energy component name set
        # {'B_41':['es', 'vdw'], ...}

        self._pair_energy_components = {}
        # {('A_32', 'B_41'):['es', 'vdw'], ('B_41', 'ligand'):[...]}
        # maps the pair energy of a position to the energy component name set

        self._energy_components = []
        # [, 'Stretch', ...], ["Advanced Electrostatics", "vdW", ...]
        # contains all energy component name set that are found in
        # any energy file

        self._total_energies = []
        # pattern: [-3213.4, -123.23, ...]
        # for each solution the total energy

        self._binding_energies = []
        # the sum of all energies of the binding scoring function
        # [-2343.0, -2346.56, ...]

        self._fold_energies = []
        # the sum of all energies of the fold energy function, i.e. not
        # involving the ligand
        # [-3453.0, -5464.0, ...]

        self.index_mapper = IndexMapper.from_index_file(index_file)
        # self.index_mapper.readIndexFile(index_file)
        self.__read(solution_file, energy_file)

    def __read(self, solution_file: str, energy_file: str) -> NoReturn:
        """
        Read the input files and fill the attributes

        Parameters
        ----------
        solution_file: str
            Path of the output file of the solver containing the solutions
        energy_file: str
            Path of the solver input file containing the energy tables
        """
        self._read_sontag_solutions(solution_file)
        self._read_sontag_energies(energy_file)
        self._calculate_global_energies()

    def _read_sontag_solutions(self, solution_file: str) -> NoReturn:
        """
        Read solution from a result file as produced by the sontag method
         (or a concatenation of these):
         one line per solution, with the indices of the solution conformers
         separated by whitespace.

        Parameters
        ----------
        solution_file: str
            Path of the solution file

        """
        with open(solution_file, 'r') as infile:
            for line in infile:
                if line.isspace():
                    continue
                sol_indices = list(map(int, line.split()))
                assert len(sol_indices) == self.index_mapper.get_pos_count(), \
                    'Number of design positions in solution file ' \
                    'does not match number in index file'
                self._raw_conf_indices.append(sol_indices)
                sol_confs = {}
                for index, conf in enumerate(sol_indices):
                    pos = self.index_mapper.get_pos_name(index)
                    sol_confs[pos] = self.index_mapper.get_restype_and_index(pos, conf)
                self._solution_confs.append(sol_confs)

    def _read_sontag_energies(self, sontag_lambdas_file: str) -> NoReturn:
        """
        Reads the energy values that are part of the solutions

        Parameters
        ----------
        sontag_lambdas_file: str
            Path of the sontag solver input file containing the energy tables ('lambdas.txt')
        """
        # holds the indices of the design position pairs in the same
        # order as they are listed in the lamdbas file
        pairs = []
        positions = self.index_mapper.get_pos_names()
        for i in range(1, self.index_mapper.get_pos_count()):
            for j in range(0, i):
                pairs.append((j, i))
        pair_count = len(pairs)
        for solution in range(self.get_solution_number()):

            self._pairwise_energies.append({})
            self._self_energies.append([])

            with open(sontag_lambdas_file, 'r') as lambda_file:
                for index, line in enumerate(lambda_file):
                    # parse the energy file
                    if line.isspace():
                        break

                    # List containing all energies in row
                    energies = line.split()
                    if index < pair_count:
                        # a pair energy has to be read, they are listed first in
                        # the lambda files
                        (position_1, position_2) = pairs[index]
                        # TODO: exclude option:
                        #  _raw_conf_indices must be updated during exclusion of amino acids in design.
                        #       Not omitting exclusion of conformers of excluded amino acid in raw_conf_indices
                        #       lead to too many taken pair indices of conformers.
                        #           result: IndexError: list index out of range
                        index_1 = self._raw_conf_indices[solution][position_1]
                        index_2 = self._raw_conf_indices[solution][position_2]

                        pos_2 = positions[position_2]

                        energy_index = self.index_mapper.get_conf_count_for_pos(pos_2) * index_1 + index_2
                        energy = -1 * float(energies[energy_index])
                        self._pairwise_energies[solution][pairs[index]] = energy
                    else:
                        # self energies of conformers
                        self_index = index - pair_count
                        conf_index = self._raw_conf_indices[solution][self_index]
                        energy = -1 * float(energies[conf_index])
                        self._self_energies[solution].append(energy)

    def _calculate_global_energies(self) -> NoReturn:
        """
        Calculate the global energy sums (total energy, binding and packing)

        """
        lig_i = self.index_mapper.get_pos_index('ligand')
        for solution in range(self.get_solution_number()):
            # calculate total energies by summing up the individual parts
            total = 0
            binding = 0
            fold = 0
            for index, self_energy in enumerate(self._self_energies[solution]):
                total += self_energy
                if index == lig_i:
                    binding += self_energy
                else:
                    fold += self_energy
            for (position_1, position_2), pair_energy in self._pairwise_energies[solution].items():
                total += pair_energy
                if position_1 == lig_i or position_2 == lig_i:
                    binding += pair_energy
                if position_1 != lig_i and position_2 != lig_i:
                    fold += pair_energy
            self._total_energies.append(total)
            self._fold_energies.append(fold)
            self._binding_energies.append(binding)

    def read_detailed_self_energies(self) -> NoReturn:
        """
        Read detailed self energies from raw energy files, according to index from self._solution_confs
        """
        for solution in self._solution_confs:
            detailed_self_energies = {}
            for position, (residue, idx) in solution.items():
                if position == 'ligand':
                    df = pd.read_csv(os.path.join(self.work_dir, 'energies', f'{self.forcefield}_{self.rotamer_path.split("/")[-1]}', f'ligand_scaffold_{self.scorer}', f'{residue}.csv'), delimiter='\t', index_col=0)
                else:
                    df = pd.read_csv(os.path.join(self.work_dir, 'energies', f'{self.forcefield}_{self.rotamer_path.split("/")[-1]}', f'sidechain_scaffold_{self.forcefield}', f'{position}_{residue}.csv'), delimiter='\t', index_col=0)
                # ['ligand_pose', 'gauss1', 'gauss2', 'repulsion', 'hydrophobic', 'hydrogen']
                # or ['vdW', 'ES']
                energy_components = list(df.columns)  # Skip ligand-pose
                if position not in self._self_energy_components.keys():
                    # {'A_171': ['vdW', 'ES'],
                    # 'A_204': ['vdW', 'ES'],
                    # 'ligand': ['gauss1', 'gauss2', 'repulsion', 'hydrophobic', 'hydrogen']}
                    self._self_energy_components[position] = energy_components
                elif energy_components not in self._energy_components:
                    # [['vdW', 'ES'], ['gauss1', 'gauss2', 'repulsion', 'hydrophobic', 'hydrogen']]
                    self._energy_components.append(energy_components)
                detailed_self_energies[position] = np.array(df)[idx, :]
            self._detailed_self_energies.append(detailed_self_energies)  # Skip ligand_2298/ASP_24 column

    def read_detailed_pair_energies(self) -> NoReturn:
        """
        Read detailed pair energies from raw energy files according to row and column indices from self._solution_confs
        """
        for solution in self._solution_confs:
            detailed_pair_energies = {}

            for pose_a, pose_b in self.index_mapper.get_pairs():
                # Switch pose variables if ligand comes first
                if pose_a == 'ligand':
                    pose_a, pose_b = pose_b, pose_a
                # Get resname corresponding to solution
                res_a, row_id = solution[pose_a]
                res_b, col_id = solution[pose_b]

                if pose_b == 'ligand':
                    df = pd.read_csv(os.path.join(self.work_dir, 'energies', f'{self.forcefield}_{self.rotamer_path.split("/")[-1]}', f'ligand_sidechain_{self.scorer}', f'{res_b}_{pose_a}_{res_a}.csv'), delimiter='\t', index_col=0)
                else:
                    df = pd.read_csv(os.path.join(self.work_dir, 'energies', f'{self.forcefield}_{self.rotamer_path.split("/")[-1]}', f'sidechain_sidechain_{self.forcefield}', f'{pose_a}_{res_a}_{pose_b}_{res_b}.csv'), delimiter='\t', index_col=0)

                energy_components = int(len(df.columns)/(int(df.columns[-1].split('_')[1]) + 1))

                self._pair_energy_components[(pose_a, pose_b)] = ['_'.join(component.split('_')[2:]) for component in df.columns[:energy_components]]
                # Reads the correct row and then reads all the energy components in the right column
                # transposed ligand array, col_id become rows, row_id become cols
                if pose_b == 'ligand':
                    detailed_pair_energies[(pose_a, pose_b)] = np.array(df)[col_id, row_id*energy_components:row_id*energy_components+energy_components]
                else:
                    detailed_pair_energies[(pose_a, pose_b)] = np.array(df)[row_id, col_id*energy_components:col_id*energy_components+energy_components]
            self._detailed_pair_energies.append(detailed_pair_energies)  # Skip ligand_2298/ASP_24 column

    def create_design_structures(self, design_path: str) -> NoReturn:
        """
        Creates structures for design solutions

        Parameters
        ----------
        design_path: Output path were to write the design structures to
        rotamer_path: Path to sidechain rotamers
        """

        # Read in ligand poses
        ligand = Molecule(self.ligand)
        ligand_poses = Molecule(self.ligand_poses)
        ligand_poses.read(f'{os.path.splitext(self.ligand_poses)[0]}.xtc')
        ligand_poses.set('segid', 'L')
        ligand_poses.set('resname', 'MOL')
        # Read in minimized protein structure
        scaffold = Molecule(self.scaffold)

        for num, solution in enumerate(self._solution_confs):

            design_structure = scaffold.copy()
            ligand_idx = solution['ligand'][1]
            lig_pose = ligand.copy()
            lig_pose.set('segid', 'L')
            lig_pose.set('resname', 'MOL')
            # Set coords to coords of ligand pose of solution
            lig_pose.set('coords', ligand_poses.coords[:, :, ligand_idx])
            for pose in self.get_positions()[:-1]:
                chain_id, res_id = pose.split('_')
                residue, rotamer = solution[pose]
                # Read in rotamers for all design positions
                residue_mol = Molecule(
                    os.path.join(self.rotamer_path, pose, f'{residue}.pdb')
                )
                residue_mol.set('chain', chain_id)
                residue_mol.set('resid', res_id)
                # Set design position to rotamer of solution
                residue_mol.dropFrames(keep=[rotamer])
                insert_atom_idx = np.where(
                    (design_structure.resid == int(res_id)) & (design_structure.chain == chain_id))[0][0]
                residue_mol.set('segid', design_structure.segid[insert_atom_idx])
                design_structure.remove(f'chain {chain_id} and resid {res_id}', _logger=False)
                design_structure.insert(residue_mol, insert_atom_idx)
            # Write ligand pose and scaffold structure to output directory
            lig_pose.write(os.path.join(design_path, str(num), 'ligand.mol2'))
            design_structure.write(os.path.join(design_path, str(num), 'receptor.pdb'))

    def get_residue_for_positions(self, solution_index: int, position: str) -> str:
        """
        Get the residue at a design pos

        Parameters
        ----------
        solution_index: index of the solution
        position: name of position (e.g. 'B_41', 'ligand')

        Returns
        -------
        name of solution residue
        """
        return self._solution_confs[solution_index][position][0]

    def get_positions(self) -> List[str]:
        """
        Get names of the design positions

        Returns
        -------
        list of names of backbone positions in design, ordered naturally
        (e.g. ['A_54', 'A_58', 'B_'41', 'ligand', ...] )
        """
        return self.index_mapper.get_pos_names()

    def get_solution_number(self) -> int:
        """

        Returns
        -------
        Number of solutions
        """
        return len(self._raw_conf_indices)

    def get_total_energy(self, solution_index: int) -> float:
        """
        Get the total energy of a solution

        Parameters
        ----------
        solution_index: index of solution

        Returns
        -------
        Total energy (float)
        """
        return self._total_energies[solution_index]

    def get_binding_energy(self, solution_index: int) -> float:
        """
        Get the total binding energy of a solution

        Parameters
        ----------
        solution_index: index of solution

        Returns
        -------
        Float of binding energy
        """
        return self._binding_energies[solution_index]

    def get_fold_energy(self, solution_index: int) -> float:
        """
        Get the total fold energy of a solution

        Parameters
        ----------
        solution_index: index of solution

        Returns
        -------
        Float of packing energy
        """
        return self._fold_energies[solution_index]

    def get_self_energy(self, solution_index: int, pos: str) -> float:
        """
        Get self energy at a position

        Parameters
        ----------
        solution_index: index of solution
        pos: backbone position of conformer (e.g. "B_41")

        Returns
        -------
        Float of energy
        """
        index = self.get_positions().index(pos)
        return self._self_energies[solution_index][index]

    def get_self_energy_component(self, solution_index: int, pos: str, component: str) -> float:
        """
        Get conformer self energy contribution of one energy component

        Parameters
        ----------
        solution_index: index of solution
        pos: design position of conformer (e.g. "B_41", "ligand")
        component: component name

        Returns
        -------
        Float of energy
        """
        index = self._self_energy_components[pos].index(component)
        return self._detailed_self_energies[solution_index][pos][index]

    def get_pair_energy(self, solution_index: int, pos1: str, pos2: str) -> float:
        """
        Get a pairwise energy

        Parameters
        ----------
        solution_index: index of solution
        pos1: first design position
        pos2: second design position

        Returns
        -------
        Float of energy
        """
        index_pos_1 = self.index_mapper.get_pos_index(pos1)
        index_pos_2 = self.index_mapper.get_pos_index(pos2)
        if (index_pos_1, index_pos_2) in self._pairwise_energies[solution_index]:
            return self._pairwise_energies[solution_index][(index_pos_1, index_pos_2)]
        return self._pairwise_energies[solution_index][(index_pos_1, index_pos_2)]

    def get_pair_energy_component(self, solution_index: int, component: str, pos1: str, pos2: str) -> float:
        """
        Get contribution of a energy component to a pairwise energy

        Parameters
        ----------
        solution_index: index of solution
        component: name of energy component
        pos1: first design position
        pos2: second design position

        Returns
        -------
        Float of energy
        """
        # if pos1 < pos2:
        #     key = (pos1, pos2)
        # else:
        key = (pos1, pos2)
        index = self._pair_energy_components[key].index(component)
        return self._detailed_pair_energies[solution_index][key][index]

    def get_mutable_position(self) -> List[str]:
        """
        Return mutable backbone positions (i.e. where more than one residue type is allowed

        Returns
        -------
        Naturally sorted list of design positions e.g ["A_57", "B_41"]
        """
        mutable_positions = []
        for pose in self.index_mapper.get_pos_names():
            if len(self.index_mapper.get_resnames(pose)) > 1:
                mutable_positions.append(pose)
        return mutable_positions

    def is_mutable(self, position: str) -> bool:
        """
        Checks if a design position is mutable, i.e. if there is more than one residue type allowed

        Parameters
        ----------
        position: backbone position name

        Returns
        -------
        True if position is mutable, False otherwise
        """
        return len(self.index_mapper.get_resnames(position)) > 1

    def get_self_energy_component_names(self, pos: str) -> List[str]:
        """
        Get names of energy components for self energies at a position

        Parameters
        ----------
        pos: design position name

        Returns
        -------
        A list of the names of the energy components of the self energy scoring function of pos
        """
        return self._self_energy_components[pos]

    def get_pair_energy_component_names(self, pos1: str, pos2: str) -> List[str]:
        """
        Get energy component names of pair energies between two positions

        Parameters
        ----------
        pos1: name of position 1
        pos2: name of position 2

        Returns
        -------
        A list of the names of energy components of the energy function used for a pairwise energy calculation
        """
        # if pos1 < pos2:
        #     return self._pair_energy_components[(pos1, pos2)]
        return self._pair_energy_components[(pos1, pos2)]

    def has_detailed_self_energies(self) -> bool:
        """
        Check if there are detailed energies (i.e. values of the individual components) available for the self energies

        Returns
        -------
        True if available, False otherwise
        """
        return self._self_energy_components != []

    def has_detailed_pair_energies(self) -> bool:
        """
        Check if there are detailed energies (i.e. values of the individual components) available for the pair energies

        Returns
        -------
        True if available, False otherwise
        """
        return self._pair_energy_components != []

    def has_detailed_self_energy(self, pos: str) -> bool:
        """
        Check if there are detailed self energies available for a position

        Parameters
        ----------
        pos: design position name

        Returns
        -------
        True if energies of the individual components are present for the design position
        """
        return pos in self._self_energy_components

    def get_pair_energy_component_sets(self) -> List[List[str]]:
        """
        Get the sets of energy component names for pair energies

        Returns
        -------
        List of lists, each list is a set of energy component names for pair energies
        """
        self._pair_energy_components.values()
