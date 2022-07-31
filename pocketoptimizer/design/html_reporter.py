import os
from typing import List, Tuple, Dict, NoReturn
import logging

import numpy as np
import pandas as pd
from natsort import natsorted
import matplotlib.pyplot as plt

import pocketoptimizer.utility.html_tags as t
from pocketoptimizer.design.design_solution import DesignSolution

logger = logging.getLogger(__name__)


class HtmlReporter:
    """
    Class for creating a set of HTML files with energy tables and plots presenting results of receptor design calculations.
    """

    def __init__(self, design_solutions: DesignSolution, sidechain_positions: Dict[str, List[str]], output_dir: str):
        """
        Constructor method.

        Parameters
        ----------
        design_solutions: DesignSolution
            The DesignSolution object with the info that should be displayed & reported
        sidechain_positions: dict
            Dictionary where key is the design position 'ChainID_ResID'
            and values are the amino acid three letter codes for this position
        output_dir: str
            The directory where the report files should be created
        """
        self._solutions = design_solutions
        self._sidechain_positions = natsorted(sidechain_positions)
        self._output_dir = output_dir

    @staticmethod
    def _get_color_val(e: float, min_e: float, max_e: float) -> str:
        """
        Determines the rgb color value of an energy value, if it is negative, it is blue, else red.
        The shade depends on its value relative to the min and max energy values.

        Parameters
        ----------
        e: float
            Energy value
        min_e: float
            Minimal energy of e's energy set
        max_e: float
            Maximal energy of e's energy set

        Returns
        -------
        String with value for HTML's style attribute to color the HTML element accordingly.
        """
        color = 'background-color:rgb('

        if e == 0 or min_e == max_e:
            return color + '255,255,255)'
        if min_e < 0 and max_e < 0:
            # adjust range
            factor = 1 - (e - max_e) / (min_e - max_e)
        elif min_e > 0 and max_e > 0:
            # adjust range
            factor = 1 - (e - min_e) / (max_e - min_e)
        else:
            factor = 1 - abs(e) / (max(abs(min_e), abs(max_e)))
        c = str(int(255. * factor))
        if e > 0:
            color += '255,' + c + ',' + c + ')'
        else:
            color += c + ',' + c + ',255)'
        return color

    def _get_table_html_tag(self, value_matrix: np.ndarray, title_rows: List[List[Tuple[str, int]]], title_columns: List[List[str]],
                            tc_highlights: List[List[int]] = None, column_groups: List[Tuple[int]] = None, add_total: bool = False) -> t.TABLE:
        """
        Returns the HTML <table> tag of the table.

        Parameters
        ----------
        value_matrix: np.ndarray
            A numpy array of the energy values
        title_rows: list
            A list of lists, here inner lists are title ROWS
            [[('Position', 1), ('Energy, '3')], [], ...].
            First Tuple  is title, second is how
            many columns should be spanned
        title_columns: list
            List of lists, inner lists are title COLUMNS
            (e.g. [['A_32:VAL', 'A_34:SER',...]]. First column
            is the one below the title row(s)
        tc_highlights: list
            List of lists, where each list corresponds to a column
            list in title_columns, and contains the indices of cells
            that should be highlighted in this column, e.g. [[1,4,..]]
        column_groups: list
            List of tuples with indexes of columns in
            value_matrix that should be grouped (and colored) together
            (i.e. (1, 3) means that colums in range(1,3) should be
            colored together, i.e. columns with index 1 and 2
        add_total: bool
            If True, a row containing total energies of the
            columns containing values is added.

        Returns
        -------
        A HTML table object (HTMLTags.TABLE)
        """
        if column_groups is None:
            column_groups = []
        column_groups.sort()
        title_column_cells = []
        for ic, col in enumerate(title_columns):
            cells = []
            for it, title in enumerate(col):
                cell = t.TD(title)
                if tc_highlights and it in tc_highlights[ic]:
                    cell.add_attributes(style='background-color:grey')
                cells.append(cell)
            title_column_cells.append(cells)
        table = t.TABLE(cellpadding='5', border='2', rules='groups')

        # group title columns
        table.insert(t.COLGROUP(span=str(len(title_column_cells))))
        for (c1, c2) in column_groups:
            # add column groups
            table.insert(t.COLGROUP(span=str(c2 - c1)))
        thead = t.THEAD()
        for row in title_rows:
            trow = t.TR()
            for (title, span) in row:
                trow.insert(t.TH(title, colspan=span, align='left'))
            thead.insert(trow)
        table.insert(thead)

        rows = [t.TR() for dummy in range(value_matrix.shape[0])]
        # insert the column containing row titles into the table:
        for col in title_column_cells:
            for i, cell in enumerate(col):
                rows[i].insert(cell)

        # get the maximal and minimal values for each column
        bounds = {}
        for (c1, c2) in column_groups:
            maxmin = (np.min(value_matrix[:, c1:c2]),
                      np.max(value_matrix[:, c1:c2]))
            for i in range(c1, c2):
                bounds[i] = maxmin

        if add_total:
            # if a line with the summed totals is requested, create it here
            trow = t.TR(t.TH('Total', colspan=len(title_column_cells)))

        for i, c in enumerate(value_matrix.transpose()):
            # loops over the matrix colums (i.e. the lines of the transpose)
            for j, v in enumerate(c):
                if i in bounds:
                    # if i is in bounds, this column shall be colored
                    if not np.isnan(v) and not np.isnan(bounds[i][0]) and not np.isnan(bounds[i][1]):
                        color = self._get_color_val(v, bounds[i][0], bounds[i][1])
                        rows[j].insert(t.TD('%.4f' % v, align='right', style=color))
                else:
                    rows[j].insert(t.TD('%.4f' % v, align='right'))
            if add_total:
                trow.insert(t.TD('%.4f' % sum(c), align='right'))
        tbody = t.TBODY()
        for row in rows:
            tbody.insert(row)
        table.insert(tbody)
        if add_total:
            table.insert(t.TBODY(trow))
        return table

    def _write_energy_overview(self, solution_index: int, body: t.BODY) -> NoReturn:
        """
        Write small table with an overview over the energies of a design solution.

        Parameters
        ----------
        solution_index: int
            Index of solution
        body: t.BODY
            Html BODY tag in which to insert the table
        """
        table = t.TABLE(cellpadding="5")
        table.insert(
            t.TR(t.TH('Total energy:', align='left'),
                 t.TD('%.2f kcal/mol' % (self._solutions.get_total_energy(solution_index)),
                      align='right')),
            t.TR(t.TH(t.A('Binding energy:', href='#binding_scores'),
                      align='left'),
                 t.TD('%.2f kcal/mol' % (self._solutions.get_binding_energy(solution_index)),
                      align='right')),
            t.TR(t.TH(t.A('Packing energy:', href='#packing_scores'),
                      align='left'),
                 t.TD('%.2f kcal/mol' % (self._solutions.get_fold_energy(solution_index)),
                      align='right')))
        body.insert(table)

    def _write_mutation_report(self, solution_index: int, body: t.BODY) -> NoReturn:
        """
        Writes a HTMl table with the amino acid types at the design positions.

        Parameters
        ----------
        solution_index: int
            Index of solution
        body: t.BODY
            Html BODY tag in which to insert the table
        """
        body.insert(t.H4('Amino acids for mutable positions'))
        table = t.TABLE(border='1', cellpadding='5', rules='groups')
        table.insert(t.THEAD(t.TR(t.TH('Position', align='left'),
                                  t.TH('Amino acid', align='left'))))
        tbody = t.TBODY()

        for p in self._sidechain_positions:
            if self._solutions.is_mutable(p):
                tbody.insert(
                    t.TR(
                        t.TD(p, align='left'),
                        t.TD(self._solutions.get_residue_for_positions(solution_index, p), align='left')
                    )
                )
        table.insert(tbody)
        body.insert(table)

    def _write_simple_ligand_pair_report(self, solution_index: int, body: t.BODY) -> NoReturn:
        """
        Writes a simple table of the ligand conformer/side chain rotamer interaction energies, i.e. no energy components, only the total energy.

        Parameters
        ----------
        solution_index: int
            Index of solution
        body: t.BODY
            Html BODY tag in which to insert the table
        """
        body.insert(t.P('Interactions between ligand '
                         'and flexible side chains'))
        title_row = [[('Position', 1), ('Total energy (scaled)', 1)]]

        value_matrix = np.zeros((len(self._solutions.get_positions() - 1), 1))
        pos_column = []
        pos_highlight = []
        for i, p in enumerate(self._sidechain_positions):
            pos_column.append(p + ':' + self._solutions.get_residue_for_positions(solution_index, p))
            e = self._solutions.get_pair_energy(solution_index, p, 'ligand')
            value_matrix[i, 0] = e
            if self._solutions.is_mutable(p):
                pos_highlight.append(i)
        i = len(self._sidechain_positions)
        for p in self._solutions.get_positions():
            if p not in self._sidechain_positions and p != 'ligand':
                pos_column.append(p + ':' + self._solutions.get_residue_for_positions(solution_index, p))
                e = self._solutions.get_pair_energy(solution_index, p, 'ligand')
                value_matrix[i, 0] = e
                if self._solutions.is_mutable(p):
                    pos_highlight.append(i)
                i += 1

        color_groups = [(0, 1)]
        result_table = self._get_table_html_tag(value_matrix,
                                                title_row, [pos_column], [pos_highlight],
                                                color_groups, add_total=True)
        body.insert(result_table)

    def _write_ligand_self_report(self, solution_index: int, body: t.BODY) -> NoReturn:
        """
        Write a report of the ligand self energies (between ligand and fixed scaffold part).

        Parameters
        ----------
        solution_index: int
            Index of solution
        body: t.BODY
            Html BODY tag in which to insert the table
        """
        body.insert(t.H4('Ligand self energies'),
                    t.P('Interactions between ligand and '
                         'fixed scaffold part'))
        scaled = self._solutions.get_self_energy(solution_index, 'ligand')

        if self._solutions.has_detailed_self_energies():
            components = self._solutions.get_self_energy_component_names('ligand')
            value_matrix = np.zeros((1, len(components) + 2))
            value_matrix[0, 0] = scaled
            summa = 0.0
            title_rows = [[('', 1), ('Energy [kcal/mol]', len(components) + 2)],
                          [('Position', 1), ('Scaled', 1), ('Sum', 1)]]
            for i, c in enumerate(components):
                title_rows[1].append((c, 1))
                e = self._solutions.get_self_energy_component(solution_index, 'ligand', c)
                summa += e
                value_matrix[0, i + 2] = e
            value_matrix[0, 1] = summa
            scaffold_table = self._get_table_html_tag(
                value_matrix,
                title_rows, [['Scaffold']], None,
                [(0, 1), (1, 2), (2, len(components) + 2)]
            )
        else:
            scaffold_table = self._get_table_html_tag(
                np.array([[scaled]]),
                [[('', 1), ('Total energy (scaled)', 1)]],
                [['Scaffold']]
            )
        body.insert(scaffold_table)

    def _get_detail_lig_sidechain_table(self, solution_index: int) -> t.TABLE:
        """
        Creates a HTML table with the detailed ligand sidechain interaction energies.

        Parameters
        ----------
        solution_index: int
            Solution index

        Returns
        -------
        An HTML table tag.
        """
        component_names = None
        pos_column = []
        pos_highlight = []
        for i, p in enumerate(self._sidechain_positions):
            if not component_names:
                # assuming all ligand-sidechain energies have same components
                component_names = self._solutions.get_pair_energy_component_names(p, 'ligand')
                title_rows = [[('', 1), ('Energy [kcal/mol]', len(component_names) + 2)],
                              [('Position', 1), ('Scaled', 1), ('Sum', 1)]]
                for name in component_names:
                    title_rows[1].append((name, 1))

                value_matrix = np.zeros((len(self._sidechain_positions),
                                         len(component_names) + 2))

            res = self._solutions.get_residue_for_positions(solution_index, p)
            pos_column.append(p + ':' + res)
            if self._solutions.is_mutable(p):
                pos_highlight.append(i)

            scaled = self._solutions.get_pair_energy(solution_index, p, 'ligand')
            value_matrix[i, 0] = scaled
            summa = 0.0
            for j, c in enumerate(component_names):
                e = self._solutions.get_pair_energy_component(solution_index, c, p, 'ligand')
                value_matrix[i, j + 2] = e
                summa += e
            value_matrix[i, 1] = summa
        color_groups = [(0, 1), (1, 2), (2, len(component_names) + 2)]

        result_table = self._get_table_html_tag(value_matrix, title_rows,
                                                [pos_column], [pos_highlight],
                                                color_groups, True)
        return result_table

    def _write_detailed_ligand_pair_report(self, solution_index: int, body: t.BODY) -> None or NoReturn:
        """
        Reports the energy values of the scoring function components between ligand conformers and side chain rotamers.

        Parameters
        ----------
        solution_index: int
            Index of solution
        body: t.BODY
            Html BODY tag in which to insert the table
        """
        body.insert(t.H4('Ligand pairwise energies'),
                    t.P('Interactions between ligand '
                        'and flexible side chains'))
        # first, amino acid side chains:
        body.insert(self._get_detail_lig_sidechain_table(solution_index))

        if self._solutions.get_positions() == natsorted(self._sidechain_positions + ['ligand']):
            return

        # now, non-sidechain positions, if they exist
        positions = self._solutions.get_positions()
        positions.remove('ligand')
        if natsorted(positions) == self._sidechain_positions:
            return

        body.insert(t.P('Interactions between ligand '
                        'and other design positions'))
        component_names = None
        pos_column = []
        pos_highlight = []
        ph = 0
        value_list = []
        for p in self._solutions.get_positions():
            if p in self._sidechain_positions or p == 'ligand':
                continue
            comp_names = self._solutions.get_pair_energy_component_names(p, 'ligand')
            if comp_names != component_names:
                if component_names:
                    # write table with old components
                    value_matrix = np.array(value_list)
                    value_list = []
                    title_rows = [[('Position', 1), ('Scaled', 1), ('Sum', 1)]]
                    for name in component_names:
                        title_rows[0].append((name, 1))
                    color_groups = [(0, 1), (1, 2), (2, len(component_names) + 2)]
                    result_table = self._get_table_html_tag(value_matrix,
                                                            title_rows, [pos_column],
                                                            [pos_highlight], color_groups,
                                                            add_total=True)
                    body.insert(result_table)
                    ph = 0
                    pos_highlight = []
                    pos_column = []
                component_names = comp_names
            res = self._solutions.get_residue_for_positions(solution_index, p)
            pos_column.append(p + ':' + res)
            if self._solutions.is_mutable(p):
                pos_highlight.append(ph)
            ph += 1
            scaled = self._solutions.get_pair_energy(solution_index, p, 'ligand')
            vl = [scaled, 0.0]
            summa = 0.0
            for c in component_names:
                e = self._solutions.get_pair_energy_component(solution_index, c, p, 'ligand')
                vl.append(e)
                summa += e
            vl[1] = summa
            value_list.append(vl)

        value_matrix = np.array(value_list)
        title_rows = [[('Position', 1), ('Scaled', 1), ('Sum', 1)]]
        for name in component_names:
            title_rows[0].append((name, 1))
        color_groups = [(0, 1), (1, 2), (2, len(component_names) + 2)]
        result_table = self._get_table_html_tag(value_matrix,
                                                title_rows, [pos_column], [pos_highlight],
                                                color_groups, add_total=True)
        body.insert(result_table)

    def _write_ligand_score_report(self, solution_index: int, body: t.BODY) -> NoReturn:
        """
        Write the report for ligand interaction scores.

        Parameters
        ----------
        solution_index: int
            Index of solution
        body: t.BODY
            Html BODY tag in which to insert the table
        """
        if self._solutions.has_detailed_pair_energies():
            self._write_detailed_ligand_pair_report(solution_index, body)
        else:
            self._write_simple_ligand_pair_report(solution_index, body)
        self._write_ligand_self_report(solution_index, body)

    def _write_detailed_self_report(self, solution_index: int, body: t.BODY) -> NoReturn:
        """
        Write a report of the self energies including the contribution of each energy component.

        Parameters
        ----------
        solution_index: int
            Index of solution
        body: t.BODY
            Html BODY tag in which to insert the table
        """
        component_names = None
        title_rows = []
        position_column = []
        position_highlight = []
        ph = 0
        value_list = []

        for pos in self._solutions.get_positions():
            if pos == 'ligand':
                continue
            comp_names = self._solutions.get_self_energy_component_names(pos)
            if comp_names != component_names:
                if component_names:
                    # write table with old components
                    value_matrix = np.array(value_list)
                    value_list = []
                    title_rows.append([('Position', 1), ('Scaled', 1),
                                       ('Sum', 1)])
                    for name in component_names:
                        title_rows[1].append((name, 1))
                    color_groups = [(0, 1), (1, 2), (2, len(component_names) + 2)]
                    result_table = self._get_table_html_tag(value_matrix,
                                                            title_rows, [position_column],
                                                            [position_highlight], color_groups, True)
                    body.insert(result_table)
                    ph = 0
                    position_column = []
                    position_highlight = []
                    title_rows = []
                else:
                    title_rows = [[('', 1), ('Energy [kcal/mol]', len(comp_names) + 1)]]
                component_names = comp_names
            if self._solutions.is_mutable(pos):
                position_highlight.append(ph)
            ph += 1
            sum_e = 0.0
            res = self._solutions.get_residue_for_positions(solution_index, pos)
            position_column.append(pos + ':' + res)
            scaled = self._solutions.get_self_energy(solution_index, pos)
            vl = [scaled, 0.0]
            for c in component_names:
                e = self._solutions.get_self_energy_component(solution_index, pos, c)
                sum_e += e
                vl.append(e)
            vl[1] = sum_e
            value_list.append(vl)

        value_matrix = np.array(value_list)
        title_rows.append([('Position', 1), ('Scaled', 1), ('Sum', 1)])
        for name in component_names:
            title_rows[len(title_rows) - 1].append((name, 1))
        color_groups = [(0, 1), (1, 2), (2, len(component_names) + 2)]
        result_table = self._get_table_html_tag(value_matrix,
                                                title_rows, [position_column], [position_highlight],
                                                color_groups, add_total=True)
        body.insert(result_table)

    def _write_simple_self_report(self, solution_index: int, body: t.BODY) -> None or NoReturn:
        """
        Write report about conformer self energies, only outputting the total energy.

        Parameters
        ----------
        solution_index: int
            Index of solution
        body: t.BODY
            Html BODY tag in which to insert the table
        """
        title_rows = [[('Position', 1), ('Energy [kcal/mol]', 1)]]
        positions = self._solutions.get_positions()
        positions.remove('ligand')
        position_column = []
        position_highlight = []
        value_matrix = np.zeros((len(positions), 1))
        for i, pos in enumerate(positions):
            e = self._solutions.get_self_energy(solution_index, pos)
            res = self._solutions.get_residue_for_positions(solution_index, pos)
            position_column.append(pos + ':' + res)
            if self._solutions.is_mutable(pos):
                position_highlight.append(i)
            value_matrix[i, 0] = e
        if np.alltrue(value_matrix == 0.):
            body.insert(t.P(t.I('All zero')))
            return
        color_groups = [(0, 1)]
        table = self._get_table_html_tag(value_matrix, title_rows,
                                         [position_column], [position_highlight],
                                         color_groups, True)
        body.insert(table)

    def _write_self_energies_report(self, solution_index: int, body: t.BODY) -> NoReturn:
        """
        Reports the self energies of the solution conformers

        Parameters
        ----------
        solution_index: index of solution
        body: html BODY tag in which to insert the table


        """
        body.insert(t.H4('Side chain self energies'),
                    t.P('Interaction with fixed scaffold part'))
        if self._solutions.has_detailed_self_energies():
            self._write_detailed_self_report(solution_index, body)
        else:
            self._write_simple_self_report(solution_index, body)

    def _write_detailed_pair_report(self, solution_index: int, body: t.BODY) -> NoReturn:
        """
        Write a detailed report for conformer pairs including the energies
           of the force field components

        Parameters
        ----------
        solution_index: index of solution
        body: html BODY tag in which to insert the table

        """
        # prepape order of rows: first , sidechain pairs, then mixed pairs,
        # then non-sidechain pairs
        position_pairs = []
        for i, p1 in enumerate(self._sidechain_positions[:-1]):
            for p2 in self._sidechain_positions[i + 1:]:
                position_pairs.append((p1, p2))

        non_aa_positions = [p for p in self._solutions.get_positions()
                            if p != 'ligand' and
                            p not in self._sidechain_positions]

        if len(non_aa_positions) > 0:
            for sp in self._sidechain_positions:
                for nap in non_aa_positions:
                    position_pairs.append((sp, nap))
            if len(non_aa_positions) > 1:
                for i, p1 in enumerate(non_aa_positions[:-1]):
                    for p2 in non_aa_positions[i + 1:]:
                        position_pairs.append((p1, p2))

        component_names = None
        position_columns = [[], []]
        position_highlights = [[], []]
        ph = 0
        value_list = []
        title_rows = None

        for (p1, p2) in position_pairs:
            res1 = self._solutions.get_residue_for_positions(solution_index, p1)
            res2 = self._solutions.get_residue_for_positions(solution_index, p2)
            comp_names = self._solutions.get_pair_energy_component_names(p1, p2)
            if comp_names != component_names:
                if component_names:
                    value_matrix = np.array(value_list)
                    value_list = []
                    title_rows.append([('Position 1', 1), ('Position 2', 1),
                                       ('Scaled', 1), ('Sum', 1)])
                    for name in component_names:
                        title_rows[len(title_rows) - 1].append((name, 1))
                    color_groups = [(0, 1), (1, 2), (2, len(component_names) + 2)]
                    result_table = self._get_table_html_tag(value_matrix,
                                                            title_rows, position_columns,
                                                            position_highlights,
                                                            color_groups, add_total=True)
                    body.insert(result_table)
                    title_rows = []
                    position_columns = [[], []]
                    position_highlights = [[], []]
                    ph = 0
                else:
                    title_rows = [[('', 1), ('', 1), ('Energy [kcal/mol]', len(comp_names) + 2)]]
                component_names = comp_names
            position_columns[0].append(p1 + ':' + res1)
            position_columns[1].append(p2 + ':' + res2)
            if self._solutions.is_mutable(p1):
                position_highlights[0].append(ph)
            if self._solutions.is_mutable(p2):
                position_highlights[1].append(ph)
            ph += 1
            vl = [self._solutions.get_pair_energy(solution_index, p1, p2), 0.0]
            for c in component_names:
                vl.append(self._solutions.get_pair_energy_component(solution_index, c, p1, p2))
            vl[1] = sum(vl[2:])
            value_list.append(vl)

        value_matrix = np.array(value_list)
        if title_rows is not None:
            title_rows.append([('Position 1', 1), ('Position 2', 1),
                               ('Scaled', 1), ('Sum', 1)])
            for name in component_names:
                title_rows[len(title_rows) - 1].append((name, 1))
            color_groups = [(0, 1), (1, 2), (2, len(component_names) + 2)]
            result_table = self._get_table_html_tag(value_matrix, title_rows,
                                                    position_columns, position_highlights,
                                                    color_groups, add_total=True)
            body.insert(result_table)

    def _write_simple_pair_report(self, solution_index: int, body: t.BODY) -> NoReturn:
        """
        Write a report of the conformer pair energies; only the total score is reported.

        Parameters
        ----------
        solution_index: index of solution
        body: html BODY tag in which to insert the table


        """
        positions = self._solutions.get_positions()
        positions.remove('ligand')
        title_rows = [[('Position 1', 1), ('Position 2', 1), ('Energy', 1)]]
        pair_count = sum([len(positions) - i for i in range(1, len(positions))])
        position_columns = [[], []]
        position_highlights = [[], []]
        ph = 0
        value_matrix = np.zeros((pair_count, 1))
        pc = 0
        for i, p1 in enumerate(positions):
            res1 = self._solutions.get_residue_for_positions(solution_index, p1)
            for p2 in positions[i + 1:]:
                res2 = self._solutions.get_residue_for_positions(solution_index, p2)
                e = self._solutions.get_pair_energy(solution_index, p1, p2)
                value_matrix[pc, 0] = e
                pc += 1
                position_columns[0].append(p1 + ':' + res1)
                position_columns[1].append(p2 + ':' + res2)
                if self._solutions.is_mutable(p1):
                    position_highlights[0].append(ph)
                if self._solutions.is_mutable(p2):
                    position_highlights[1].append(ph)
                ph += 1
        color_groups = [(0, 1)]
        table = self._get_table_html_tag(value_matrix, title_rows,
                                         position_columns, position_highlights,
                                         color_groups, True)
        body.insert(table)

    def _write_pairwise_energies_report(self, solution_index: int, body: t.BODY) -> NoReturn:
        """
        Writes the pairwise energies of the solution conformer pairs

        Parameters
        ----------
        solution_index: index of solution
        body: html BODY tag in which to insert the table


        """
        body.insert(t.H4('Side chain pairwise energies'),
                    t.P('Interactions between flexible side chains'))
        if self._solutions.has_detailed_pair_energies():
            self._write_detailed_pair_report(solution_index, body)
        else:
            self._write_simple_pair_report(solution_index, body)

    def create_reports(self) -> NoReturn:
        """
        Creates a report for each solution listing the mutations and the energy values.

        """
        for s in range(self._solutions.get_solution_number()):
            head = t.HEAD(t.TITLE('Report for design solution %s' % s))
            body = t.BODY()

            navtable = t.TABLE()
            row = t.TR()
            if s != 0:
                row.insert(t.TD(t.A('&laquo; Previous solution', align='left',
                                    href='../' + str(s - 1) + '/report.html')))
            row.insert(t.TD())
            row.insert(t.TD(t.A('Overview', align='center',
                                href='../summary.html')))
            row.insert(t.TD())
            if s < self._solutions.get_solution_number() - 1:
                row.insert(t.TD(t.A('Next solution &raquo;', align='right',
                                    href='../' + str(s + 1) + '/report.html')))
            navtable.insert(row)
            body.insert(navtable)
            body.insert(t.BR())

            body.insert(t.H1('Energy report for solution %s' % s))
            self._write_energy_overview(s, body)
            body.insert(t.BR())
            self._write_mutation_report(s, body)
            body.insert(t.BR())
            body.insert(t.A(t.H2('Ligand binding scores'),
                            name='binding_scores'))
            self._write_ligand_score_report(s, body)
            body.insert(t.BR())
            body.insert(t.A(t.H2('Scaffold packing scores'),
                            name='packing_scores'))
            self._write_self_energies_report(s, body)
            body.insert(t.BR())
            self._write_pairwise_energies_report(s, body)

            outdir = os.path.join(self._output_dir, str(s))
            if not os.path.exists(outdir):
                os.makedirs(outdir, exist_ok=True)
            with open(os.path.join(outdir, 'report.html'), 'w') as out:
                print(t.HTML(head, body), file=out)

    def create_summary_plot(self) -> NoReturn:
        """
        Create a png of a bar plot with the total binding and
        total packing scores of all _solutions depicted.
        """
        binding = []
        packing = []
        for i in range(self._solutions.get_solution_number()):
            binding.append(self._solutions.get_binding_energy(i))
            packing.append(self._solutions.get_fold_energy(i))

        ind = np.arange(self._solutions.get_solution_number())

        # One plot with splited bars indicating share of binding and packing energy
        df = pd.DataFrame(index=ind, data={'Packing': packing, 'Binding': binding})
        ax = df.plot(kind='bar', stacked=True, edgecolor='black', colormap='Paired')
        # Draw black line for total energy
        df.sum(axis=1).plot(ax=ax, color='k', label='Total', linestyle="", marker="o")

        # format the x axis
        ax.set_xlabel('Solution')
        # label the y-axis (Total energy)
        ax.set_ylabel('Total Energy [kcal/mol]')

        ax.set_xlim(-1, self._solutions.get_solution_number())
        # Rotate x-axis if more than 25 solutions to make easier readable
        if len(self._solutions._raw_conf_indices) >= 25:
            ax.tick_params(axis='x', rotation=45)

        # show a legend
        ax.legend(bbox_to_anchor=(1.25, 1), loc='upper right', ncol=1)

        # save a png of the plot
        plt.savefig(os.path.join(self._output_dir, 'summary.png'), dpi=300, bbox_inches='tight')

    def create_sequence_logo(self) -> NoReturn:
        """
        Create a sequence logo of the mutations in the solutions
        """

        import weblogo as wll
        from IPython.display import Image, display
        from pocketoptimizer.preparation.aacodes import three2one

        sequences = wll.seq.SeqList()
        for si in range(self._solutions.get_solution_number()):
            seq = ''
            for p in self._solutions.get_mutable_position():
                res = self._solutions.get_residue_for_positions(si, p)
                if res in three2one:
                    seq += three2one[res]
            sequences.append(seq)
        sequences.alphabet = wll.seq.unambiguous_protein_alphabet
        data = wll.LogoData.from_seqs(sequences)
        options = wll.LogoOptions()
        options.title = 'Titulo'
        options.resolution = 350.
        options.fineprint = ''
        options.creator_text = 'Receptor Design Pipeline'
        options.color_scheme = wll.chemistry
        options.unit_name = 'probability'
        options.yaxis_label = 'Frequency'
        options.annotate = self._solutions.get_mutable_position()
        options.fontsize = 5
        options.number_fontsize = 4
        logoformat = wll.LogoFormat(data, options)
        with open(os.path.join(self._output_dir, 'seqlogo.png'), 'wb') as fout:
            fout.write(wll.png_formatter(data, logoformat))

        display(Image(os.path.join(self._output_dir, 'seqlogo.png')))

    def create_summary(self, plot: bool = True, logo: bool = True) -> NoReturn:
        """
        Create a summary html page with info of the design mutations and energies

        Parameters
        ----------
        plot: if True, a plot of the energes will be created
        logo: if True, a sequence logo for mutable positions will be created

        """
        mutable_positions = self._solutions.get_mutable_position()

        head = t.HEAD(t.TITLE('Summary'))
        body = t.BODY(t.H1('Summary of Design Solutions'),
                      t.BR())
        table = t.TABLE(t.COLGROUP(span='1'), t.COLGROUP(span='3'),
                        t.COLGROUP(span='3'),
                        border='2', cellpadding='5', rules='groups')
        headrow1 = t.TR(t.TH('Solution'),
                        t.TH('Residues', colspan=str(len(mutable_positions))),
                        t.TH('Energies [kcal/mol]', colspan='3'))
        headrow2 = t.TR(t.TH('&nbsp'))
        for p in mutable_positions:
            headrow2.insert(t.TH(p))
        headrow2.insert(t.TH('Total'), t.TH('Binding'), t.TH('Packing'))
        table.insert(t.THEAD(headrow1, headrow2))
        tbody = t.TBODY()
        for s in range(self._solutions.get_solution_number()):
            row = t.TR()
            if os.path.exists(os.path.join(self._output_dir, str(s), 'report.html')):
                row.insert(t.TD(t.A(str(s), href=str(s) + '/report.html')))
            else:
                row.insert(t.TD(str(s)))
            for p in mutable_positions:
                row.insert(t.TD(self._solutions.get_residue_for_positions(s, p)))
            row.insert(t.TD('%.2f' % self._solutions.get_total_energy(s)),
                       t.TD('%.2f' % self._solutions.get_binding_energy(s)),
                       t.TD('%.2f' % self._solutions.get_fold_energy(s)))
            tbody.insert(row)
        table.insert(tbody)
        body.insert(table)

        if plot:
            body.insert(t.BR())
            body.insert(t.H2('Total Energy Plot:'))
            logger.info('Create energy plot.')
            self.create_summary_plot()
            body.insert(t.IMG(src='./summary.png', style='max-width: 500px; max-height: 300px'))

        if logo:
            try:
                logger.info('Create sequence logo for mutable positions.')
                self.create_sequence_logo()
                body.insert(t.BR())
                body.insert(t.H2('Sequence Logo for mutable design positions:'))
                body.insert(t.BR())
                body.insert(t.IMG(src='./seqlogo.png', style='max-width: 450px; max-height: 600px'))
            except:
                logger.error('Sequence logo creation failed. Omitting sequence logo creation.')

        with open(os.path.join(self._output_dir, 'summary.html'), 'w') as out:
            print(t.HTML(head, body), file=out)


