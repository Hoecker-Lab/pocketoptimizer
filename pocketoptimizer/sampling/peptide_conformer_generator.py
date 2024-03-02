import os
import shutil
import numpy as np
import tempfile as tf
from tqdm.auto import tqdm
import multiprocessing as mp
from functools import partial
from typing import List, Dict, NoReturn
import logging
from moleculekit.molecule import Molecule
from ffevaluation.ffevaluate import FFEvaluate

from pocketoptimizer.sampling.sidechain_rotamers_ffev import FFRotamerSampler

logger = logging.getLogger(__name__)


class PeptideSampler(FFRotamerSampler):
    """ Class for peptide conformer sampling based on rotamers using either
    backbone independent rotamers from cmlib or backbone dependent rotamers from dunbrack rotamer library.
    cmlib:
        * superimpose rotamers onto structure
    dunbrack:
        * read rotamers for respective phi/psi angle combination
        * get all chi angles for aa
        * convert to radians
        * set dihedral angles
        * prune based on probability
    - save conformers into pdb file
    """
    def __init__(self, library: str, positions: List[Dict[str, str]], **kwargs):
        """Constructor method

        Parameters
        ---------
        library: str
            Rotamer library to use

        """
        super().__init__(**kwargs)
        self.library = library
        self.positions = positions

    def merge_conformers(self, conf_ids: List[int], outfile: str) -> NoReturn:
        """
        Writes out rotamers as models into one merged .pdb file

        Parameters
        ----------
        conf_ids: int
            Conformers with an energy below the threshold
        outfile: name/path of the output file
        """
        with open(outfile, 'w') as merged_rot_file:
            for conf_id in conf_ids:
                with open(os.path.join(self.tmp_dir, f'ligand_conf_{conf_id}.pdb'), 'r') as conf_file:
                    for line in conf_file:
                        if line.startswith('MODEL'):
                            line = f'MODEL        {str(conf_id)}\n'
                        # skip 'END' lines (distinguish between 'END' and 'ENDMDL')
                        elif line.startswith('END') and not line.startswith('ENDMDL'):
                            continue
                        merged_rot_file.write(line)
            merged_rot_file.write('END')

    def calculate_energy(self, conf_id: int, structure: Molecule, ffev: FFEvaluate) -> float:
        """
        Calculates the vdW energy of a peptide conformation

        Parameters
        ----------
        conf_id: int
            Conformer id to calculate
        structure: class:moleculekit.Molecule
            Object containing conformers
        ffev: :class:ffevaluation.ffevaluate.FFevaluate
            Force field object

        Returns
        -------
        Returns energy
        """
        # Set coordinates to coordinates of conformer
        structure.set('coords', structure.coords[:, :, conf_id])
        energies = ffev.calculateEnergies(structure.coords[:, :, conf_id])
        structure.write(os.path.join(self.tmp_dir, f'ligand_conf_{conf_id}.pdb'))

        return energies['vdw']

    def conformer_sampling(self, vdw_filter_thresh: float = 100.0, include_native: bool = True) -> NoReturn:
        """
        Parameters
        ----------
        vdw_filter_thresh: float
            Filtering threshold to avoid internal clashes [default: 100 kcal/mol]
        include_native: bool
            Include the starting conformation
        """
        from pocketoptimizer.utility.utils import MutationProcessor, load_ff_parameters, calculate_chunks

        outfile = os.path.join(self.work_dir, 'ligand', self.forcefield, 'conformers', 'ligand_confs.pdb')

        if os.path.isfile(outfile):
            logger.info('Conformers are already sampled.')
            return

        self.tmp_dir = tf.mkdtemp(dir=self.tmp_dir, prefix='calculateRotamers_')
        os.chdir(self.tmp_dir)

        logger.info('Start conformer sampling procedure.')

        mutation_processor = MutationProcessor(structure=self.ligand_protonated,
                                               mutations=self.positions,
                                               forcefield=self.forcefield)
        termini_positions = mutation_processor.check_termini()
        struc, prm = load_ff_parameters(structure_path=self.build_scaffold_params,
                                        forcefield=self.forcefield)
        # Remove the protein
        struc.filter('segid L', _logger=False)

        # Set coordinates of minimized peptide
        min_struc = Molecule(self.ligand_protonated)
        struc.coords = min_struc.coords

        # Create object to hold all conformers
        confs = struc.copy()

        for position in self.positions:

            chain = position['chain']
            resid = position['resid']
            resname = struc.get('resname', f'chain {chain} and resid {resid} and name CA')[0]

            # Keep original rotamer
            residue = struc.copy()
            residue.filter(f'chain {chain} and resid {resid}', _logger=False)
            if resname != 'GLY' and resname != 'ALA':
                if self.library == 'cmlib':
                    rotamers = self.read_db(resname=resname)

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
                    for i, torsion in enumerate(self.sidechain_torsions[resname]):
                        # select the four atoms forming the dihedral angle according to their atom names
                        current_rot.setDihedral([int(current_rot.get('index', sel=f'name {torsion[0]}')),
                                                 int(current_rot.get('index', sel=f'name {torsion[1]}')),
                                                 int(current_rot.get('index', sel=f'name {torsion[2]}')),
                                                 int(current_rot.get('index', sel=f'name {torsion[3]}'))],
                                                rotamer[i] * (np.pi/180), bonds=bonds)
                    # append rotameric states as frames to residue
                    residue.appendFrames(current_rot)

                if not include_native:
                    residue.dropFrames(drop=0)

            num_confs = confs.coords.shape[-1]
            for conf_id in range(num_confs):
                modified_conf = struc.copy()
                modified_conf.set('coords', confs.coords[:, :, conf_id])
                for rotamer_id in range(residue.coords.shape[-1]):
                    modified_conf.set('coords', residue.coords[:, :, rotamer_id], sel=f'chain {chain} and resid {resid}')
                    confs.appendFrames(modified_conf)
            # Remove the starting conformation
            confs.dropFrames(drop=num_confs - 1)

        nconfs = confs.coords.shape[-1]

        # Generate FFEvaluate object
        ffev = FFEvaluate(struc, prm)

        energies = np.ndarray(nconfs)
        with tqdm(total=nconfs, desc='Filter Conformer(s)') as pbar:
            with mp.Pool(processes=self.ncpus) as pool:
                for conf_id, energy in enumerate(pool.imap(
                        partial(self.calculate_energy,
                                structure=confs.copy(),
                                ffev=ffev
                                ), np.arange(nconfs),
                        chunksize=calculate_chunks(nposes=nconfs, ncpus=self.ncpus))):
                    energies[conf_id] = energy
                    pbar.update()
        val_ids = [val_id[0] for val_id in np.argwhere(energies <= min(energies) + vdw_filter_thresh)]

        self.merge_conformers(conf_ids=val_ids, outfile=outfile)
        logger.info(f'Generated {len(val_ids)} conformer(s).')

        if os.path.isdir(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)