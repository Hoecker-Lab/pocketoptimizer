import os
from typing import List, Dict, NoReturn
from natsort import natsorted
from io import StringIO

from pocketoptimizer.design.design_solution import DesignSolution


class TxtReporter:
    """
    Class for creating a set of text files containing energy tables and mutation results of a set of design solutions.
    """

    def __init__(self, design_solutions: DesignSolution, sidechain_positions: Dict[str, List[str]], output_dir: str):
        """
        Constructor method.

        Parameters
        ----------
        design_solutions: DesignSolution
            Object of the design solutions that should be displayed & reported
        sidechain_positions: Dictionary
            Key is the design position 'ChainID_ResID'
            and values are the amino acid three letter codes for this position
        output_dir: str
            Output directory where the report files should be created
        """
        self._sidechain_positions = natsorted(sidechain_positions)
        self._solutions = design_solutions
        self._output_dir = output_dir

    def _write_energy_overview(self, solution_index: int, out_file: StringIO) -> NoReturn:
        """
        Write a short overview of the solution energy.

        Parameters
        ----------
        solution_index: int
            Index of solution
        out_file: StringIO
            File object of the output file
        """
        out_file.write(' TOTAL ENERGY:    ' +
                       ('%.2f kcal/mol' % (self._solutions.get_total_energy(solution_index)))
                       .rjust(10) + '\n BINDING ENERGY:  ' +
                       ('%.2f kcal/mol' % (self._solutions.get_binding_energy(solution_index)))
                       .rjust(10) + '\n PACKING ENERGY:  ' +
                       ('%.2f kcal/mol' % (self._solutions.get_fold_energy(solution_index)))
                       .rjust(10) + '\n\n')

    def _write_mutation_report(self, solution_index: int, out_file: StringIO) -> NoReturn:
        """
        Writes a report with the amino acid types at the design positions.

        Parameters
        ----------
        solution_index: int
            Index of solution
        out_file: StringIO
            File object of the output file
        """

        out_file.write(' AMINO ACIDS:\n\n'
                       ' for all mutable positions\n\n'
                       '   Position | Amino acid\n   ---------------------\n')

        for p in self._sidechain_positions:
            if self._solutions.is_mutable(p):
                out_file.write('    ' + p.ljust(13) +
                               self._solutions.get_residue_for_positions(solution_index, p)
                               + '\n')
        out_file.write('\n')

    def _write_detailed_ligand_pair_report(self, solution_index: int, out_file: StringIO) -> NoReturn:
        """
        Writes reports for energies between ligand conformer and sidechain rotamer, including all energy components of the respective scoring function.

        Parameters
        ----------
        solution_index: int
            Index of solution
        out_file: StringIO
            File object of the output file
        """

        col_width = 0
        for p in self._solutions.get_positions():
            res = self._solutions.get_residue_for_positions(solution_index, p)
            col_width = max(col_width, len(res) + len(p) + 1)

        out_file.write('  ' + 'Position'.ljust(col_width) + ' | Energies [kcal/mol]\n')
        title_line = '  ' + ' ' * col_width + ' |  Scaled ||     Sum'

        totals = {}
        totals['sum'] = 0.0
        totals['scaled'] = 0.0
        # assuming all ligand-sidechain energies have same components
        comp_names = self._solutions.get_self_energy_component_names('ligand')
        for c in comp_names:
            totals[c] = 0.0
            title_line += ' | ' + c.ljust(8)
        out_file.write(title_line + '\n  ' + '-' * len(title_line) + '\n')

        # first, the amino acid side chains:
        for p in self._sidechain_positions:
            res = self._solutions.get_residue_for_positions(solution_index, p)
            scaled = self._solutions.get_pair_energy(solution_index, p, 'ligand')
            totals['scaled'] += scaled
            summa = 0.0
            out_str = ''
            for c in comp_names:
                e = self._solutions.get_pair_energy_component(solution_index, c, p, 'ligand')
                totals[c] += e
                summa += e
                out_str += ('%.4f' % e).rjust(max(8, len(c))) + '   '
            out_file.write('  ' + (p + ':' + res).ljust(col_width) + '   ' +
                           ('%.4f' % scaled).rjust(8) + '   ' +
                           ('%.4f' % summa).rjust(8) + '   ' + out_str + '\n')
            totals['sum'] += summa

        out_file.write('\n   ' + 'Total'.rjust(col_width) + '  '
                       + ('%.4f' % totals['scaled']).rjust(8) +
                       '   ' + ('%.4f' % totals['sum']).rjust(8))
        for c in comp_names:
            out_file.write('   '
                           + ('%.4f' % totals[c]).rjust(max(8, len(c))))
        out_file.write('\n\n')

    def _write_simple_ligand_pair_report(self, solution_index: int, out_file: StringIO) -> NoReturn:
        """
        Write a simple report of the ligand conformer/side chain rotamer energies, no energy components, only the total energy.

        Parameters
        ----------
        solution_index: int
            Index of solution
        out_file: StringIO
            File object of the output file
        """
        out_file.write('             |  Total energy (scaled)'
                       '   ----------------------------------\n')
        for p in self._sidechain_positions:
            scaled = self._solutions.get_pair_energy(solution_index, p, 'ligand')
            out_file.write('  ' + (p + ':' + self._solutions.get_residue_for_positions(solution_index, p)).rjust(10)
                           + '   '
                           + ('%.4f' % scaled).rjust(11) + '\n')
        for p in self._solutions.get_positions():
            if p not in self._sidechain_positions and p != 'ligand':
                scaled = self._solutions.get_pair_energy(solution_index, p, 'ligand')
                out_file.write('  ' + (p + ':' + self._solutions.get_residue_for_positions(solution_index, p)).rjust(10)
                               + '   ' + ('%.4f' % scaled).rjust(11) + '\n')
        out_file.write('\n')

    def _write_ligand_self_report(self, solution_index: int, out_file: StringIO) -> NoReturn:
        """
        Write a report of the ligand self energies (i.e. ligand-fixed scaffold).

        Parameters
        ----------
        solution_index: int
            Index of solution
        out_file: StringIO
            File object of the output file
        """
        if self._solutions.has_detailed_self_energies():
            comp_names = self._solutions.get_self_energy_component_names('ligand')
            out_file.write('\n')
            title_row = '             |  Scaled ||     Sum'
            for c in comp_names:
                title_row += ' | ' + c.rjust(8)
            out_file.write(title_row + '\n   ' + '-' * len(title_row) + '\n')

            scaled = self._solutions.get_self_energy(solution_index, 'ligand')
            summa = 0.0
            out_str = ''
            for c in self._solutions.get_self_energy_component_names('ligand'):
                e = self._solutions.get_self_energy_component(solution_index, 'ligand', c)
                summa += e
                out_str += ('%.4f' % e).rjust(max(8, len(c))) + '   '
            out_file.write('             ' + ('%.4f' % scaled).rjust(10) + ' ' +
                           ('%.4f' % summa).rjust(10) + '   ' + out_str + '\n')
        else:
            out_file.write('   Total energy:   ' +
                           ('%.4f' % self._solutions.get_self_energy(solution_index, 'ligand')).rjust(10) + '\n')

    def _write_ligand_score_report(self, solution_index: int, out_file: StringIO) -> NoReturn:
        """
        Write the report for ligand interaction scores.

        Parameters
        ----------
        solution_index: int
            Index of solution
        out_file: StringIO
            File object of the output file
        """

        out_file.write('LIGAND INTERACTION SCORES \n\n')
        out_file.write('SELF ENERGIES \n')
        out_file.write('(Interactions between ligand and fixed scaffold part)\n\n')
        self._write_ligand_self_report(solution_index, out_file)
        out_file.write('\nPAIRWISE ENERGIES \n')
        out_file.write('(Interactions between ligand and the '
                       'flexible sidechains)\n\n')
        if self._solutions.has_detailed_pair_energies():
            self._write_detailed_ligand_pair_report(solution_index, out_file)
        else:
            self._write_simple_ligand_pair_report(solution_index, out_file)
        out_file.write('\n')

    def _write_detailed_self_report(self, solution_index: int, out_file: StringIO) -> NoReturn:
        """
        Write a report of the self energies including the contribution of each energy component.

        Parameters
        ----------
        solution_index: int
            Index of solution
        out_file: StringIO
            File object of the output file
        """
        col_width = 0
        for pos in self._solutions.get_positions():
            if pos != 'ligand':
                comp_names = self._solutions.get_self_energy_component_names(pos)
                res = self._solutions.get_residue_for_positions(solution_index, pos)
                col_width = max(col_width, len(res) + len(pos) + 1)

        totals = {}
        totals['sum'] = 0.0
        totals['scaled'] = 0.0
        out_file.write('  ' + 'Position'.ljust(col_width) + ' | Energies [kcal/mol]')
        title_row = '  ' + ' '.rjust(col_width) \
                    + ' |   Scaled  ||       Sum'

        for c in comp_names:
            totals[c] = 0.0
            title_row += ' | ' + c.rjust(8)
        title_row += '\n'
        out_file.write('\n' + title_row + '  ' + '-' * len(title_row) + '\n')
        for pos in self._solutions.get_positions():
            if pos != 'ligand':
                res = self._solutions.get_residue_for_positions(solution_index, pos)
                out_file.write('   ' + (pos + ':' + res).ljust(col_width)
                               + '  ')
                scaled = self._solutions.get_self_energy(solution_index, pos)
                totals['scaled'] += scaled
                summa = 0.0
                out_str = ''
                for c in comp_names:
                    e = self._solutions.get_self_energy_component(solution_index, pos, c)
                    summa += e
                    totals[c] += e
                    out_str += ('%.4f' % e).rjust(max(8, len(c))) + '   '
                totals['sum'] += summa
                out_file.write(('%.4f' % scaled).rjust(10) + '   '
                               + ('%.4f' % summa).rjust(10)
                               + '   ' + out_str + '\n')

        out_file.write('\n      ' + 'Total'.rjust(col_width) + '  '
                       + ('%.4f' % totals['scaled']).rjust(8) +
                       '   ' + ('%.4f' % totals['sum']).rjust(8))
        for c in comp_names:
            out_file.write('   '
                           + ('%.4f' % totals[c]).rjust(max(8, len(c))))
        out_file.write('\n\n')

    def _write_simple_self_report(self, solution_index: int, out_file: StringIO) -> NoReturn:
        """
        Write report about conformer self energies, only outputting the total energy

        Parameters
        ----------
        solution_index: int
            Index of solution
        out_file: StringIO
            File object of the output file.
        """
        out_file.write('   Position  | Energy [kcal/mol]')
        out_file.write('   -------------------------\n')
        for pos in self._solutions.get_positions():
            if pos == 'ligand':
                continue
            e = self._solutions.get_self_energy(solution_index, pos)
            res = self._solutions.get_residue_for_positions(solution_index, pos)
            out_file.write('   ' + (pos + ':' + res).ljust(11))
            out_file.write(('%.4f' % e).rjust(10) + '\n')
        out_file.write('\n')

    def _write_detailed_pair_report(self, solution_index: int, out_file: StringIO) -> NoReturn:
        """
        Write a detailed report for conformer pairs including the energies of the force field components.

        Parameters
        ----------
        solution_index: int
            Index of solution.
        out_file: StringIO
            File object of the output file.
        """

        # prepape order of rows: first, sidechain pairs, then mixed pairs,
        # then non-sidechain pairs
        position_pairs = []
        for i, p1 in enumerate(self._sidechain_positions[:-1]):
            for p2 in self._sidechain_positions[i + 1:]:
                position_pairs.append((p1, p2))

        non_aa_positions = self._solutions.get_positions()
        non_aa_positions.remove('ligand')
        for sp in self._sidechain_positions:
            non_aa_positions.remove(sp)

        if len(non_aa_positions) > 0:
            for sp in self._sidechain_positions:
                for nap in non_aa_positions:
                    position_pairs.append((sp, nap))
            if len(non_aa_positions) > 1:
                for i, p1 in enumerate(non_aa_positions[:-1]):
                    for p2 in non_aa_positions[i + 1:]:
                        position_pairs.append((p1, p2))

        col_width = 0
        for p in self._solutions.get_positions():
            if p != 'ligand':
                comp_names = self._solutions.get_self_energy_component_names(p)
                res = self._solutions.get_residue_for_positions(solution_index, p)
                col_width = max(col_width, len(res) + len(p) + 1)

        totals = {}
        totals['sum'] = 0.0
        totals['scaled'] = 0.0
        out_file.write('  ' + 'Position 1'.ljust(col_width) + ' | ' +
                       'Position 2'.ljust(col_width) + ' | Energies [kcal/mol]\n')
        title_row = '  ' + ' ' * col_width + ' | ' + ' ' * col_width + \
                    ' |   Scaled ||      Sum  '
        for c in comp_names:
            totals[c] = 0.0
            title_row += ' | ' + c.rjust(8)
        out_file.write(title_row + '\n  ' + '-' * len(title_row) + '\n')

        for (p1, p2) in position_pairs:
            res1 = self._solutions.get_residue_for_positions(solution_index, p1)
            res2 = self._solutions.get_residue_for_positions(solution_index, p2)

            out_file.write('  ' + (p1 + ':' + res1).ljust(col_width) + '   ' + (p2 + ':' + res2).ljust(col_width))
            scaled = self._solutions.get_pair_energy(solution_index, p1, p2)
            totals['scaled'] += scaled
            summa = 0.0
            out_str = ''
            for c in comp_names:
                e = self._solutions.get_pair_energy_component(solution_index, c, p1, p2)
                summa += e
                totals[c] += e
                out_str += ('%.4f' % e).rjust(max(8, len(c)) + 2) + ' '
            totals['sum'] += summa
            out_file.write(('%.4f' % scaled).rjust(10) + '   ' +
                           ('%.4f' % summa).rjust(10) + '   ' + out_str + '\n')

        out_file.write('\n             ' + 'Total'.rjust(col_width) + '  '
                       + ('%.4f' % totals['scaled']).rjust(8) +
                       '   ' + ('%.4f' % totals['sum']).rjust(8))
        for c in comp_names:
            out_file.write('   '
                           + ('%.4f' % totals[c]).rjust(max(8, len(c))))
        out_file.write('\n')

    def _write_simple_pair_report(self, solution_index: int, out_file: StringIO) -> NoReturn:
        """
        Write a report of the conformer pair energies; only the total score is reported.

        Parameters
        ----------
        solution_index: int
            Index of solution
        out_file: StringIO
            File object of the output file

        """

        out_file.write('   Position 1 | Position 2 | Energy [kcal/mol]\n'
                       '   --------------------------------\n')
        lines = []
        positions = self._solutions.get_positions()
        positions.remove('ligand')
        for i, p1 in enumerate(positions):
            res1 = self._solutions.get_residue_for_positions(solution_index, p1)
            for p2 in positions[i + 1:]:
                res2 = self._solutions.get_residue_for_positions(solution_index, p2)
                line = '   ' + (p1 + ':' + res1).ljust(13) + (p2 + ':' + res2).ljust(12)
                e = self._solutions.get_pair_energy(solution_index, p1, p2)
                # _pairwise_energies pattern: [{(3, 5): 10.23, ...},...]
                line += ('%.4f' % e).rjust(10) + '\n'
                lines.append(line)
        lines.sort()
        out_file.writelines(lines)
        out_file.write('\n')

    def _write_protein_score_report(self, solution_index: int, out_file: StringIO) -> NoReturn:
        """
        Write the report for protein interaction scores.

        Parameters
        ----------
        solution_index: int
            Index of solution
        out_file: StringIO
            File object of the output file.
        """
        out_file.write('PROTEIN INTERACTION SCORES \n\n')
        out_file.write(' SELF ENERGIES\n (Interactions '
                       'between flexible sidechains and fixed scaffold part)\n\n')
        if self._solutions.has_detailed_self_energies():
            self._write_detailed_self_report(solution_index, out_file)
        else:
            self._write_simple_self_report(solution_index, out_file)

        out_file.write(' PAIRWISE ENERGIES\n (Interactions '
                       'between flexible sidechains and fixed scaffold part)'
                       '\n\n')
        if self._solutions.has_detailed_pair_energies():
            self._write_detailed_pair_report(solution_index, out_file)
        else:
            self._write_simple_pair_report(solution_index, out_file)

    def create_reports(self) -> NoReturn:
        """
        Creates a report for each solution listing the mutations and the energy values.
        """
        for s in range(self._solutions.get_solution_number()):
            outdir = os.path.join(self._output_dir, str(s))
            if not os.path.exists(outdir):
                os.makedirs(outdir, exist_ok=True)

            with open(os.path.join(outdir, 'report.txt'), 'w') as out:
                out.write('\n')
                self._write_energy_overview(s, out)
                out.write('_' * 100 + '\n\n')
                self._write_mutation_report(s, out)
                out.write('_' * 100 + '\n\n')
                self._write_ligand_score_report(s, out)
                out.write('_' * 100 + '\n\n')
                self._write_protein_score_report(s, out)

    def create_summary(self) -> NoReturn:
        """
        Create a summary text file with info of the design mutations and energies
        """
        mutable_positions = self._solutions.get_mutable_position()

        with open(os.path.join(self._output_dir, 'summary.txt'), 'w') as out:
            out.write('\n  SUMMARY OF DESIGN SOLUTIONS\n\n')
            # column width only works if position name lengths
            # are >= residue name lengths
            col_width = sum(map(len, mutable_positions)) + 3 * len(mutable_positions) - 3
            out.write('  SOLUTION | ' + 'RESIDUES'.ljust(col_width)
                      + ' | ENERGIES [kcal/mol] \n')
            out.write('           |')
            for p in mutable_positions:  # attern: [{'B_41':('ARG', 23),...}]
                out.write(' ' + p + ' |')
            out.write('   Total   |   Binding   |  Packing \n')
            out.write('  ' + '-' * (col_width + 48) + '\n')
            for s in range(self._solutions.get_solution_number()):
                out.write(str(s).rjust(8) + '     ')
                for p in mutable_positions:
                    out.write(self._solutions.get_residue_for_positions(s, p).ljust(len(p)) + '   ')
                out.write(('%.2f' % self._solutions.get_total_energy(s)).ljust(10)
                          + '   '
                          + ('%.2f' % self._solutions.get_binding_energy(s)).ljust(10)
                          + '   '
                          + ('%.2f' % self._solutions.get_fold_energy(s)).ljust(10)
                          + '\n')
            out.close()
