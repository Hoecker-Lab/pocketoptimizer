import os
from typing import NoReturn
from io import StringIO

from pocketoptimizer.utility.index_mapper import IndexMapper
from pocketoptimizer.utility.energy_reader import EnergyReader


class SontagWriter:
    """
    Class for creating the necessary input files for the sontag solver.

    Created on Nov 16, 2010

    @author: malisi


    based on analysis of "save_gmplp_state_for_c.m" output.
    this is probably only valid if all residue/pose and residue/residue pairs for
    all positions interact. if this is not true for a problem, further analysis of
    the matlab analysis has to be done, or noras "ours2yanover.m" etc.
    has to be used.
    "save_gmplp_state_for_c.m" does the conversion from matlab to c input files.
    It creates 5 files:
    * regions.txt:
     * contains one line for each position pair,
       the position indices starting at 1, not 0
     * sorted by the second index as priority (e.g. 2 3 is before 1 4, because 3<4)
     * at the end, contains a line for each position with its index
    * lambdas.txt:
     * contains the energies
     * first [number of position pairs] lines contain the pairwise energies,
       lines are ordered as in regions.txt
     * the last [number of positions] lines contain the self energies
     * energy values in a line separated by spaces
     * all energy values negated
     * energies in a row seem to be ordered 'normally'
    * intersects.txt:
     * same as regions.txt
    * var_sizes.txt
     * contains one line per position,
       containing only the number of rotamers for this position
    * region_intersects.txt
     * number of lines is number of position pairs + number of positions
     * first column counts up from 1 in all lines
     * for the first [number of position pairs] lines, second and third column
       behave as the columns in regions.txt, but not starting from 1,
       instead starting from [number of position pairs +1]
     * the last [number of positions] lines have only the first column
    """

    def __init__(self, index_mapper: IndexMapper, energy_reader: EnergyReader):
        """
        Constructor method.

        Parameters
        ----------
        index_mapper: :class: pocketoptimizer.utility.index_mapper.IndexMapper object
            Containing the design positions, residues, conformer names and indices
        energy_reader: :class: pocketoptimizer.utility.energy_reader.EnergyReader object
            Containing the energy values of the designs
        """
        self.index_mapper = index_mapper
        self.energies = energy_reader
        self.pair_index_order = []

    def write_sontag(self, output_dir: str) -> NoReturn:
        """
        Writes the files in the format necessary for the sontag solver.

        Parameters
        ----------
        output_dir: str
            Path of the directory where the files are to be written
        """
        self._write_regions(output_dir)
        self._write_varsizes(output_dir)
        self._write_regions_intersects(output_dir)
        with open(os.path.join(output_dir, 'lambdas.txt'), 'w') as lambda_file:
            self._write_pairwise_energies(lambda_file)
            self._write_self_energies(lambda_file)

    def _write_regions(self, output_dir: str) -> NoReturn:
        """
        Writes the 'regions.txt' and 'intersects.txt' files.

        Parameters
        ----------
        output_dir: str
            Path to where files should be created
        """
        with open(os.path.join(output_dir, 'regions.txt'), 'w') as regions:
            with open(os.path.join(output_dir, 'intersects.txt'), 'w') as intersects:
                # +1 because of indices have to start at 1 in the output
                for i in range(2, self.index_mapper.get_pos_count() + 1):
                    for j in range(1, i):
                        regions.write(f'{j} {i} \n')
                        intersects.write(f'{j} {i} \n')
                for i in range(1, self.index_mapper.get_pos_count() + 1):
                    regions.write(f'{i} \n')
                    intersects.write(f'{i} \n')

    def _write_varsizes(self, output_dir: str) -> NoReturn:
        """
        Writes file 'var_sizes.txt'.

        Parameters
        ----------
        output_dir: str
            Path to where files should be created
        """
        with open(os.path.join(output_dir, 'var_sizes.txt'), 'w') as f:
            for pos in self.index_mapper.get_pos_names():
                f.write(str(self.index_mapper.get_conf_count_for_pos(pos)) + '\n')

    def _write_regions_intersects(self, output_dir: str) -> NoReturn:
        """
        Writes the file 'region_intersects'.

        Parameters
        ----------
        output_dir: str
            Path to where files should be created
        """
        with open(os.path.join(output_dir, 'region_intersects.txt'), 'w') as f:
            pos_count = self.index_mapper.get_pos_count()
            pair_count = int(((pos_count * pos_count) - pos_count) / 2)
            i = 1
            for j in range(pair_count + 2, pos_count + pair_count + 1):
                for k in range(pair_count + 1, j):
                    f.write(f'{i} {k} {j}\n')
                    i += 1
            for j in range(i, i + pos_count):
                f.write(f'{j} \n')

    def _write_self_energies(self, sontag_file: StringIO) -> NoReturn:
        """
        Writes the sontag 'lambdas.txt' input section that contains the
        self energies of the poses and rotamers (i.e. the energies with
        the fixed scaffold part).

        Parameters
        ----------
        sontag_file: :class: io.StringIO object
            Writable file object for the lambda.txt file
        """
        for p in self.index_mapper.get_pos_names():
            for ic in range(self.index_mapper.get_conf_count_for_pos(p)):
                e = -1 * self.energies.get_self_energy(p, ic)
                sontag_file.write(f'{e} ')
            sontag_file.write('\n')
            sontag_file.flush()

    def _write_pairwise_energies(self, sontag_file: StringIO) -> NoReturn:
        """
        Writes the sontag 'lambda.txt' input section that contains the
        pairwise energies of poses and side chain rotamers.

        Parameters
        ----------
        sontag_file: :class: io.StringIO object
            Writable file object for the lambda.txt file
        """
        positions = self.index_mapper.get_pos_names()
        for i, pos_i in enumerate(positions):
            for j in range(0, i):
                pos_j = positions[j]
                for cj in range(self.index_mapper.get_conf_count_for_pos(pos_j)):
                    for ci in range(self.index_mapper.get_conf_count_for_pos(pos_i)):
                        e = -1 * self.energies.get_pair_energy(pos_j, pos_i, cj, ci)
                        sontag_file.write(f'{e} ')
                sontag_file.write('\n')
                sontag_file.flush()

