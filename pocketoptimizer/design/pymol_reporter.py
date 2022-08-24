import os
from typing import List, Dict, Union, NoReturn
from natsort import natsorted

from pocketoptimizer.design.design_solution import DesignSolution


class PymolReporter:
    """
    Provides methods to create pymol scripts for design solutions that
    visualize the designed structures by e.g. highlighting design positions.
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
        if not self._solutions.peptide:
            self._lig_ftype = 'mol2'
        else:
            self._lig_ftype = 'pdb'

    def create_solution_script(self, solution_index: int, scaffold_pdb: str, wt_ligand: str = None, waters: Dict[str, str] = None,
                               metals: List[List[Union[str, float]]] = None, cofactors=None) -> NoReturn:
        """
        Creates a pymol script for a design solution that loads the pdb file
           of the design and the wt, selects and displays the design positions,
           colored by their contribution to the binding energy.

        Parameters
        ----------
        solution_index: int
            The index of the solution (from the original
            kingsford solution log, is the same in the DesignSolution object)
        scaffold_pdb: str
            Path of the PDB file of the WT scaffold
        wt_ligand: str
            The path the structure file of the wt ligand. if
            None, no WT ligand representation is included in script [default: None]
        waters: dict
            Dictionary of explicit waters in design, keys are the design
            position names of the waters, values the paths to the corresponding
            PDBs [default: None]
        metals: list
            List of lists [[elem, x,y,z],[..],...] [default: None]
        cofactors: [default: None]
        """
        outdir = os.path.join(os.path.abspath(self._output_dir), str(solution_index))
        pocket = ''
        mutable = ''
        for i, pos in enumerate(self._sidechain_positions):
            assert pos in self._solutions.get_positions(), f'Position: {pos} is no valid design position.'
            chain, index = pos.split('_')
            if i > 0:
                pocket += ' or '
            pocket += '(chain ' + chain + ' and resi ' + index + ' and not' \
                      + ' (name C or name N or name O or name H or name HA))'
            if self._solutions.is_mutable(pos):
                if mutable != '':
                    mutable += ' or '
                mutable += '(chain ' + chain + ' and resi ' + index + ')'

        # Adjust colors of ligand elements, because of custom atomtypes assigned with Parameterize
        commands = [
            'bg_color white',
            'load ' + outdir + os.sep + 'receptor.pdb, design',
            'load ' + outdir + os.sep + f'ligand.{self._lig_ftype}, ligand',
            'load ' + os.path.abspath(scaffold_pdb) + ', WT_scaffold',
            'hide lines',
            'cartoon loop',
            'color gray80, WT_scaffold and elem C',
            'color gray90, design and elem C',
            'set stick_radius, 0.15',
            'show stick, ligand',
            'hide cartoon ligand',
            'select old_pocket, (' + pocket + ') and WT_scaffold',
            'select new_pocket, (' + pocket + ') and design',
            'select mutable, (' + mutable + ') and design',
            'label mutable and name ca and design, "%s - %s" %(resn,resi)',
            'show stick, old_pocket',
            'show lines, WT_scaffold and not old_pocket',
            'hide lines, WT_scaffold and name c+o+n',
            'show stick, new_pocket',
            'set_bond stick_radius, 0.25, mutable',
            'spectrum b, red_blue, new_pocket',
            'color bluewhite, ligand and elem C',
            'color splitpea, WT_scaffold and elem C',
            'color gray80, WT_scaffold and not old_pocket and elem C',
            'select don, (elem N or elem O and (neighbor hydro))',
            'select acc, (elem O or (elem N and not (neighbor hydro)))']
        if waters:
            s_create = 'create wt_waters, '
            s_delete = 'delete '
            for i, w in enumerate(waters):
                commands.append('load %s, water%i' % (waters[w], i))
                s_create += 'water%i or ' % i
                s_delete += 'water%i or ' % i
            commands.append(s_create[:-4])
            commands.append(s_delete[:-4])
            commands += [
                'show sphere, wt_waters',
                'color gray80, wt_waters',
                'select waters, resn HOH and design',
                'show sphere, waters',
                'color red, waters',
                'set sphere_scale, 0.2']

        if metals:
            metal_str = ''
            for m in metals:
                metal_str += ' %s' % m[0]  # m[0} should be the Element name
            commands += [
                'select metals, resn %s and design' % metal_str,
                'show sphere, metals',
                'set sphere_scale, 0.5']
        if cofactors:
            cf_str = ''
            for cf, pfile in cofactors.items():
                commands += ['load %s, %s' % (pfile, cf)]
                cf_str += '%s or ' % cf
            cf_str = cf_str[:-4]
            commands += ['show sticks, %s' % cf_str]
            commands += ['color yellow, (%s) and elem C' % cf_str]
            commands += ['hide everything, %s and elem H' % cf_str]

        #    bondables = []
        #    if waters:
        #        bondables.append('waters')
        #    if metals:
        #        bondables.append('metals')
        #    if cofactors:
        #        for cf in cofactors:
        #            bondables.append(cf)
        # b_str = ' or '.join(bondables)
        # TODO: create proper polar interactions

        if waters and metals:
            commands += [
                'dist HBA, ((waters or metals or ligand) and acc), '
                '((waters or metals or design) and don), 3.2',
                'dist HBD, ((waters or metals or ligand) and don), '
                '(waters or metals or design) and acc), 3.2']
        elif metals:
            commands += [
                'dist HBA, ((metals or ligand) and acc), (design and don), 3.2',
                'dist HBD, ((metals or ligand) and don), (design and acc), 3.2']
        elif waters:
            commands += [
                'dist HBA, ((waters or ligand) and acc), '
                '((waters or design) and don), 3.2',
                'dist HBD, ((waters or ligand) and don), '
                '((waters or design) and acc), 3.2']
        else:
            commands += [
                'dist HBA, (ligand and acc), (design and don), 3.2',
                'dist HBD, (ligand and don), (design and acc), 3.2']
        commands += [
            'delete don',
            'delete acc',
            'hide labels, HBA',
            'hide labels, HBD',
            'color lightorange, HBA',
            'color lightorange, HBD',
            'set label_size, 20']
        if wt_ligand:
            commands += [
                'load ' + os.path.abspath(wt_ligand) + ', WT_ligand',
                'show stick, WT_ligand',
                'color white, WT_ligand and elem C']
        commands += [
            'zoom new_pocket',
            'hide (hydro)',
            'deselect']
        with open(os.path.join(outdir, 'design.pml'), 'w') as pml:
            pml.write('\n'.join(commands) + '\n')

    def create_summary_script(self, scaffold_pdb: str, wt_ligand: str = None, waters: Dict[str, str] = None):
        #                        metals: List[List[Union[str, float]]] = None, cofactors = None):
        """
        Creates a pymol script for all design solutions that loads the pdb files
        of the designs and the wt, selects and displays the design positions.

        Parameters
        ----------
        scaffold_pdb: str
            Path of the PDB file of the WT scaffold
        wt_ligand: str
            The path the structure file of the wt ligand.
            If None, no WT ligand representation is included in script [default: None]
        waters: dict
            Dictionary of explicit waters in design, keys are the design
            position names of the waters, values the paths to the corresponding PDBs
        """
        outdir = os.path.abspath(self._output_dir)
        pocket = ''
        mutable = ''

        for i, pos in enumerate(self._sidechain_positions):
            assert pos in self._solutions.get_positions(), f'Position: {pos} is no valid design position.'
            chain, index = pos.split('_')
            if i > 0:
                pocket += ' or '
            pocket += '(chain ' + chain + ' and resi ' + index + ' and not' \
                      + ' (name C or name N or name O or name H or name HA))'
            if self._solutions.is_mutable(pos):
                if mutable != '':
                    mutable += ' or '
                mutable += '(chain ' + chain + ' and resi ' + index + ')'
        commands = ['load ' + os.path.abspath(scaffold_pdb) + ', WT_scaffold',
                    'load %s, WT_ligand' % wt_ligand,
                    'show stick, WT_ligand',
                    'hide cartoon, WT_ligand']

        for sol_i in range(self._solutions.get_solution_number()):
            dir_i = os.path.join(outdir, str(sol_i))
            commands += [
                f'load %s%sligand.{self._lig_ftype}, ligand_poses' % (dir_i, os.sep),
                'load %s%sreceptor.pdb, designs' % (dir_i, os.sep)]
        if waters:
            s_create = 'create wt_waters, '
            s_delete = 'delete '
            for i, w in enumerate(waters):
                commands += [
                    'load %s, water%i' % (waters[w], i)
                ]
                s_create += 'water%i or ' % i
                s_delete += 'water%i or ' % i
            commands.append(s_create[:-4])
            commands.append(s_delete[:-4])
            commands += [
                'show sphere, wt_waters',
                'color gray80, wt_waters',
                'select waters, resn HOH and design',
                'show sphere, waters',
                'color red, waters',
                'set sphere_scale, 0.2']

        commands += [
            'bg_color white',
            'hide lines',
            'cartoon loop',
            'set stick_radius, 0.15',
            'show stick, ligand_poses',
            'hide cartoon, ligand_poses',
            'select old_pocket, (' + pocket + ') and WT_scaffold',
            'select new_pocket, (' + pocket + ') and designs',
            'select mutable, (' + mutable + ') and designs',
            'label mutable and name ca and design, "%s - %s" %(resn,resi)',
            'show stick, old_pocket',
            'show stick, new_pocket',
            'show lines, WT_scaffold and not old_pocket',
            'hide lines, WT_scaffold and name c+o+n',
            'color gray90, designs and elem C',
            'color gray80, WT_scaffold and elem C',
            'color bluewhite, new_pocket and elem C',
            'set_bond stick_radius, 0.25, mutable',
            'set label_size, 20']
        commands += ['select don, (elem N or elem O and (neighbor hydro))',
                     'select acc, (elem O or (elem N and not (neighbor hydro)))',
                     'color bluewhite, ligand_poses and elem C',
                     'color gray80, WT_ligand and elem C'
                     ]
        if waters:
            commands += [
                'dist HBA, ((waters or ligand) and acc), '
                '((waters or design) and don), 3.2',
                'dist HBD, ((waters or ligand) and don), '
                '((waters or design) and acc), 3.2']
        else:
            commands += [
                'dist HBA, (ligand_poses and acc), (designs and don), 3.2',
                'dist HBD, (ligand_poses and don), (designs and acc), 3.2']
        commands += [
            'delete don',
            'delete acc',
            'hide labels, HBA',
            'hide labels, HBD',
            'color lightorange, HBA',
            'color lightorange, HBD']

        if waters:
            commands += [
                'show sphere, wt_waters',
                'color gray80, wt_waters',
                'select waters, resn HOH and design',
                'show sphere, waters',
                'color red, waters',
                'set sphere_scale, 0.2']
        # TODO: add cofactors
        with open(os.path.join(outdir, 'summary.pml'), 'w') as pml:
            pml.write('\n'.join(commands) + '\n')
            pml.write('zoom new_pocket\nhide (hydro)\ndeselect\n')

    def create_pymol_scripts(self, scaffold_pdb: str, wt_ligand: str = None, waters: Dict[str, str] = None,
                             metals: List[List[Union[str, float]]] = None, cofactors=None):
        """
        Create pymol scripts for visualizing the design binding pockets

        Parameters
        ----------
        scaffold_pdb: str
            Path to the PDB file of the WT scaffold
        wt_ligand: str
            A structure file of the WT ligand, if none, no representation of WT ligand is included in scripts [default: None]
        waters: dict
            Dictionary of explicit waters in design, keys are the design
            position names of the waters, values the paths to the corresponding
            PDBs [default: None]
        metals: list
            List of lists [[elem, x,y,z],[..],...] [default: None]
        cofactors: [default: None]

        """
        for sol_i in range(self._solutions.get_solution_number()):
            self.create_solution_script(sol_i, scaffold_pdb, wt_ligand, waters,
                                   metals, cofactors)
        self.create_summary_script(scaffold_pdb, wt_ligand, waters)
        # metals, cofactors)
