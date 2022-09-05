import os
import numpy as np
from tqdm.auto import tqdm
import multiprocessing as mp
from functools import partial
import tempfile as tf
import shutil
from typing import List, Tuple, Dict, Union, NoReturn
import logging

from moleculekit.molecule import Molecule
from ffevaluation.ffevaluate import FFEvaluate

logger = logging.getLogger(__name__)


class FFRotamerSampler:
    """ Class for force field based rotamer sampling using either backbone independent rotamers
    from cmlib or backbone dependent rotamers from dunbrack rotamer library.
        cmlib:
            * superimpose rotamers onto structure
            * filter based on energy
        dunbrack:
            * read rotamers for respective phi/psi angle combination
            * get all chi angles for aa
            * convert to radians
            * set dihedral angles
            * prune based on probability
            * filter based on energy
    - save rotamers into pdb files and energies into csv files
    """
    def __init__(self, work_dir: str, mutations: List[Dict[str, Union[str, List[str]]]], forcefield: str,
                 rot_path: str, library: str = 'dunbrack', tmp: str = '/tmp/'):
        """
        Constructor method.

        Parameters
        ----------
        work_dir: str
            Path to working directory
        mutations: list
            List of all mutations where rotamers should be sampled for. Contains dictionaries as list entries.
            The dictionary should contain 'chain', 'resid', 'mutations' as keys.
        forcefield: str
            Force field to use for energy calculations
        rot_path: str
            Output path for sampled rotamers
        library: str
            Library to use for selecting rotamers, either setting coordinates from pdb after superimposing (cmlib)
            or setting dihedral angles from dunbrack [default: 'dunbrack']
        tmp: str
            Temporary directory to write temporary files [default: '/tmp']
        """
        self.work_dir = work_dir
        self.mutations = mutations
        self.library = library
        self.forcefield = forcefield
        self.rot_path = rot_path
        self.tmp = tmp

    def read_db(self, resname: str, phi_angle: float = 0.0, psi_angle: float = 0.0, prob_cutoff: float = -1) -> Dict[str, List[Tuple[float]]]:
        """
        Reading function for rotamer libraries.

        Parameters
        ----------

        resname: str
            Name of the residue table to read
        phi_angle: float
            Backbone phi angle of residue in degree [default: 0.0]
        psi_angle: float
            Backbone psi angle of residue in degree [default: 0.0]
        prob_cutoff: float
            Rotamers occouring with probability less than prob_cutoff will
            be pruned if their rotameric mode is not unique [default: -1] (no pruning)

        Returns
        -------
            Dictionary for rotamers containing chi angles and standard deviations
            for rotamers belonging to a respective phi/psi angle combination if available
        """
        import sqlite3 as sl
        import pocketoptimizer.path as po_path

        def roundup(x: float, n: int = 10) -> int:
            """
            Rounds a float number to nearest n integer number.

            Parameters
            ----------
            x: float
                Number to round
            n: int
                Integer number to roundup [default: 10]

            Returns
            -------
            Nearest n integer number.
            If nearest integer number is 360, than return 0,
            since the backbone dihedral angles are periodically recurring (360° = 0°)
            """
            import math
            res = math.ceil(x/n)*n
            if (x % n < n/2) and (x % n > 0):
                res -= n
            if res == 360:
                return 0
            else:
                return res

        rotamers = {'chi': [],
                    'std': []}
        # Connect to rotamer library
        if self.library == 'dunbrack':
            con_dunbrack = sl.connect(os.path.join(po_path.path(rotamers=True), 'dunbrack.db'))
            with con_dunbrack:
                data = con_dunbrack.execute(f'SELECT * FROM {resname} WHERE phi = {roundup(phi_angle)} AND psi = {roundup(psi_angle)} AND prob >= {prob_cutoff}')
            for rotamer in data:
                rotamers['chi'].append(tuple(rotamer[4:8]))
                rotamers['std'].append(tuple(rotamer[8:]))
        elif self.library == 'cmlib':
            con_cmlib = sl.connect(os.path.join(po_path.path(rotamers=True), 'cmlib.db'))
            with con_cmlib:
                data = con_cmlib.execute(f'SELECT * FROM {resname}')
            for rotamer in data:
                rotamers['chi'].append(tuple(rotamer[1:]))
        else:
            logger.error(f'Could not read from library: {self.library}.')
            raise ValueError(f'Could not read from library: {self.library}.')

        return rotamers

    @staticmethod
    def expand_dunbrack(rotamers: List[Dict[str, float or List[int or float]]], expand: List[str] = ['chi1', 'chi2']) -> List[Tuple[float]]:
        """
        Expands defined chi-angles by +/- 2 Std

        Parameters
        ----------
        rotamers: list
            List of dictionaries for rotamers containing probability, rotameric mode, chi angles and standard deviations
            for rotamers belonging to a respective phi/psi angle combination
        expand: list
            List of which chi angles to expand, [default: ['chi1', 'chi2']]

        Returns
        -------
        List of Chi-angle tuples
        """

        # List of lists for each rotamer containing a tuple for
        # every chi-angle of an amino acid with its possible values
        rotamer_chi_angles = []

        for i, rotamer in enumerate(rotamers['chi']):
            rotamer_chi_angles.append([])
            for j, chi_angle in enumerate(rotamer):
                if f'chi{str(j+1)}' in expand:
                    rotamer_chi_angles[i].append([chi_angle - 2 * rotamers['std'][i][j],
                                                  chi_angle,
                                                  chi_angle + 2 * rotamers['std'][i][j]])
                else:
                    rotamer_chi_angles[i].append([chi_angle])

        rotamers = {'chi': []}

        for possible_rotamer in rotamer_chi_angles:
            for chi_1 in possible_rotamer[0]:
                for chi_2 in possible_rotamer[1]:
                    for chi_3 in possible_rotamer[2]:
                        for chi_4 in possible_rotamer[3]:
                            rotamer = (chi_1, chi_2, chi_3, chi_4)
                            if not rotamer in rotamers['chi']:
                                rotamers['chi'].append(rotamer)
        return rotamers

    def calculate_vdw(self, rot_id: int, structure: Molecule, ffev: FFEvaluate,
                      res_coords: np.ndarray, resname: str, resid: str, chain: str) -> np.float:
        """
        Calculates vdW energy of a rotamer

        Parameters
        ----------
        rot_id: int
            Rotamer ID
        structure: class:moleculekit.Molecule
            Object containing rotamers
        ffev: :class:ffevaluation.ffevaluate.FFevaluate
            Object between rotamer and scaffold.
        res_coords: np.array
            Rotamer coordinates.
        resname: str
            residue name
        resid: str
            Residue ID.
        chain: str
            Chain ID.

        Returns
        -------
        Returns vdw energy
        """
        # Set coordinates of side chain to coordinates from rotamer
        structure.set('coords', res_coords[:, :, rot_id], f'resid {resid} and chain {chain}')
        energies = ffev.calculateEnergies(structure.coords)
        structure.write(os.path.join(self.tmp, f'{chain}_{resid}_{resname}_rot_{rot_id}.pdb'),
                        sel=f'resid {resid} and chain {chain}')

        return energies['vdw']

    def merge_rotamers(self, rot_ids: List[int], resname: str, resid: str, chain: str, outfile: str) -> NoReturn:
        """
        Writes out rotamers as models into one merged .pdb file

        Parameters
        ----------
        rot_ids: list
            List of rotamer_ID files to merge
        resname: str
            Residue name
        resid: str
            Residue ID.
        chain: str
            Chain ID.
        outfile: name/path of the output file
        """
        with open(outfile, 'w') as merged_rot_file:
            i = 0
            for rot_id in rot_ids:
                with open(os.path.join(self.tmp, f'{chain}_{resid}_{resname}_rot_{rot_id}.pdb'), 'r') as rot_file:
                    for line in rot_file:
                        if line.startswith('MODEL'):
                            i += 1
                            line = f'MODEL        {str(i)}\n'
                        elif line.startswith('ATOM'):
                            line = line[:22] + str(rot_id).rjust(4) + line[26:]
                        elif line.startswith('CONECT'):
                            continue
                        # skip 'END' lines (distinguish between 'END' and 'ENDMDL')
                        elif line.startswith('END') and not line.startswith('ENDMDL'):
                            continue
                        merged_rot_file.write(line)
            merged_rot_file.write('END')

    def rotamer_sampling(self, vdw_filter_thresh: float = 100.0, dunbrack_filter_thresh: float = -1,
                         expand: List[str] = ['chi1', 'chi2'], ncpus: int = 1, _keep_tmp: bool = False) -> NoReturn:
        """
        Parameters
        ----------
        vdw_filter_thresh: float
            Filtering threshold to avoid clashes. [default: 100 kcal/mol]
        dunbrack_filter_thresh: float
            Filter threshold, rotamers having probability of occurence lower than filter threshold will
            be pruned if their rotameric mode does occur more than once
            (-1: no pruning, 1: pruning of all rotamers with duplicate rotamer modes) [default: -1]
        expand: list
            List of chi angles to expand [default: ['chi1', 'chi2']]
        ncpus: int
            Number of CPUs to use for multiprocessing. [default: 1]
        _keep_tmp: bool
            If the tmp directory should be deleted or not. Useful for debugging. [default: False]
        """
        from pocketoptimizer.utility.utils import MutationProcessor, load_ff_parameters, write_energies, calculate_chunks
        from pocketoptimizer.utility.molecule_types import backbone_atoms

        logger.info('Start rotamer sampling procedure using FFEvaluate.')
        if self.forcefield == 'amber_ff14SB':
            from pocketoptimizer.utility.molecule_types import _SIDECHAIN_TORSIONS_AMBER as _SIDECHAIN_TORSIONS
        elif self.forcefield == 'charmm36':
            from pocketoptimizer.utility.molecule_types import _SIDECHAIN_TORSIONS_CHARMM as _SIDECHAIN_TORSIONS
        else:
            logger.error('Force field not implemented.')
            raise NotImplementedError('Force field not implemented.')

        self.tmp = tf.mkdtemp(dir=self.tmp, prefix='calculateRotamers_')
        os.chdir(self.tmp)
        logger.info(f"Using {ncpus} CPU's for multiprocessing.")

        mutation_processor = MutationProcessor(structure=os.path.join(self.work_dir, 'scaffold', self.forcefield, 'scaffold.pdb'),
                                               mutations=self.mutations,
                                               forcefield=self.forcefield)
        termini_positions = mutation_processor.check_termini()

        for mutation in self.mutations:
            resid = mutation['resid']
            chain = mutation['chain']

            for resname in mutation['mutations']:
                # Check computation status of rotamer
                outfile = os.path.join(self.rot_path,  f'{chain}_{resid}', f'{resname}.pdb')
                if os.path.exists(outfile):
                    logger.info(f'Rotamers for residue: {chain}_{resid}_{resname} already sampled.')
                    continue
                logger.info(f'Rotamers for residue: {chain}_{resid}_{resname} not sampled yet.')

                structure_path = os.path.join(self.work_dir, 'scaffold', self.forcefield,
                                              'protein_params', f'{chain}_{resid}_{resname}')

                # Check if directories containing single mutated scaffolds are existing
                if not os.path.isfile(os.path.join(structure_path, 'structure.pdb')):
                    logger.error(f'Missing single mutated structure for mutation: {chain}_{resid}_{resname}.')
                    raise FileNotFoundError(f'Missing single mutated structure for mutation: {chain}_{resid}_{resname}.')

                struc, prm = load_ff_parameters(structure_path=structure_path, forcefield=self.forcefield)

                # Keep original rotamer
                residue = struc.copy()
                residue.filter(f'chain {chain} and resid {resid}', _logger=False)
                if resname != 'GLY' and resname != 'ALA':
                    if self.library == 'cmlib':
                        rotamers = self.read_db(resname=resname)
                    elif self.library == 'dunbrack':
                        # Check if the position is at the N-or C-terminus of the protein
                        if f'{chain}_{resid}' in termini_positions:
                                self.library = 'cmlib'
                                rotamers = self.read_db(resname=resname)
                                self.library = 'dunbrack'
                        else:
                            phi_angle = struc.getDihedral([int(struc.get('index', sel=f'chain {chain} and resid {str(int(resid)-1)} and name C')),
                                                           int(struc.get('index', sel=f'chain {chain} and resid {resid} and name N')),
                                                           int(struc.get('index', sel=f'chain {chain} and resid {resid} and name CA')),
                                                           int(struc.get('index', sel=f'chain {chain} and resid {resid} and name C'))
                                                           ]) * (180/np.pi) + 180

                            psi_angle = struc.getDihedral([int(struc.get('index', sel=f'chain {chain} and resid {resid} and name N')),
                                                           int(struc.get('index', sel=f'chain {chain} and resid {resid} and name CA')),
                                                           int(struc.get('index', sel=f'chain {chain} and resid {resid} and name C')),
                                                           int(struc.get('index', sel=f'chain {chain} and resid {str(int(resid)+1)} and name N'))
                                                           ]) * (180/np.pi) + 180

                            # Read histidine rotamers for different HIS protonation states
                            if resname in ['HID', 'HIE', 'HIP', 'HSD', 'HSE', 'HSP']:
                                _resname = 'HIS'
                            else:
                                _resname = resname

                            rotamers = self.read_db(resname=_resname,
                                                    phi_angle=phi_angle,
                                                    psi_angle=psi_angle,
                                                    prob_cutoff=dunbrack_filter_thresh)

                            rotamers = self.expand_dunbrack(rotamers=rotamers,
                                                            expand=expand)

                    else:
                        logger.error(f'Library: {self.library} not implemented.')
                        raise ValueError(f'Library: {self.library} not implemented.')

                    # Keep original rotamer
                    current_rot = residue.copy()
                    # Break proline rings, since no dihedral setting possible
                    if resname == 'PRO':
                        current_rot.deleteBonds(sel='name N or name CD', inter=False)
                    bonds = current_rot.bonds
                    # Iterate over all rotamers
                    for rotamer in rotamers['chi']:
                        for i, torsion in enumerate(_SIDECHAIN_TORSIONS[resname]):
                            # select the four atoms forming the dihedral angle according to their atom names
                            current_rot.setDihedral([int(current_rot.get('index', sel=f'name {torsion[0]}')),
                                                     int(current_rot.get('index', sel=f'name {torsion[1]}')),
                                                     int(current_rot.get('index', sel=f'name {torsion[2]}')),
                                                     int(current_rot.get('index', sel=f'name {torsion[3]}'))],
                                                    rotamer[i] * (np.pi/180), bonds=bonds)
                        # append rotameric states as frames to residue
                        residue.appendFrames(current_rot)

                nrots = residue.coords.shape[-1]

                fixed_selection = f'not segid L and not (chain {chain} and resid {resid})'
                variable_selection = f'chain {chain} and resid {resid} and not ({backbone_atoms})'
                # Generate FFEvaluate object for scoring between sidechain and rest of protein without sidechain and backbone from sidechain
                ffev = FFEvaluate(struc, prm, betweensets=(fixed_selection, variable_selection))

                energies = np.ndarray(nrots)
                with tqdm(total=nrots, desc=f'{chain}_{resid}_{resname}') as pbar:
                    with mp.Pool(processes=ncpus) as pool:
                        for i, energy in enumerate(pool.imap(
                                partial(self.calculate_vdw,
                                        structure=struc.copy(),
                                        ffev=ffev,
                                        res_coords=residue.coords,
                                        resname=resname,
                                        resid=resid,
                                        chain=chain), np.arange(nrots),
                                chunksize=calculate_chunks(nposes=nrots, ncpus=ncpus))):
                            energies[i] = energy
                            pbar.update()

                val_ids = [val_id[0] for val_id in np.argwhere(energies <= min(energies) + vdw_filter_thresh)]
                logger.info(f'Writing {len(val_ids)}/{nrots} rotamers within energy threshold of {str(vdw_filter_thresh)} kcal/mol for {resname} at position: {chain}_{resid}.')

                os.makedirs(os.path.split(outfile)[0], exist_ok=True)
                rotamer_energy_file = os.path.join(os.path.split(outfile)[0], f'{resname}.csv')

                write_energies(outpath=rotamer_energy_file,
                               energies=energies,
                               energy_terms=['Vdw'],
                               name_a=resname,
                               nconfs_a=nrots)

                self.merge_rotamers(rot_ids=val_ids, resname=resname, resid=resid, chain=chain, outfile=outfile)

        if not _keep_tmp:
            if os.path.isdir(self.tmp):
                shutil.rmtree(self.tmp)
        logger.info('Rotamer sampling procedure is finished.')
