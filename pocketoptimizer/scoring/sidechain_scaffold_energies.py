import multiprocessing as mp
import os
from functools import partial
from typing import Tuple, NoReturn
import logging
import numpy as np
from ffevaluation.ffevaluate import FFEvaluate
from moleculekit.molecule import Molecule
from tqdm.auto import tqdm

from pocketoptimizer.utility.utils import Storer

logger = logging.getLogger(__name__)


class SidechainSelfScorer(Storer):
    """
    Scorer class for sidechain/scaffold interactions
    Utilizing molecular dynamics binding energy related computations
    Requires built complexes
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.energy_terms = ['Vdw', 'Elec']

    def calculate_energy(self, id: int, struc: Molecule, ffev: FFEvaluate,
                         res_coords: np.ndarray, resid: str, chain: str) -> np.ndarray:
        """
        Calculating the energy between a rotamer and a fixed protein scaffold, or
        calculate the internal energy

        Parameters
        ----------
        id: int
            Index of current rotamer
        struc: :class: moleculekit.molecule.Molecule
            Object containing all rotamers as frames
        ffev: :class: ffevaluation.ffevaluate.FFEvaluate
            Object created for respective forcefield from parameter file
            between sets fixed_selection and variable_selection,
            where energies can be calculated from
        res_coords: np.ndarray
            Coordinates of the rotamers
        resid: str
            Residue ID of the mutated residue
        chain: str
            Chain ID of the mutated residue

        Returns
        -------
        Array containing vdw and electrostatic energies
        """

        # Set coordinates to current rotamer
        struc.set('coords', res_coords[:, :, id], f'chain {chain} and resid {resid}')
        energies = ffev.calculateEnergies(struc.coords)
        return np.array([energies['vdw'], energies['elec'] * 0.01])

    def calculate_scaffold(self) -> NoReturn:
        """
        Wrapper function to calculate sidechain/scaffold energies using
        molecular mechanics force field energy calculations
        """
        from pocketoptimizer.utility.utils import load_ff_parameters, write_energies, calculate_chunks
        from pocketoptimizer.utility.molecule_types import backbone_atoms

        logger.info(f'Compute energies using forcefield: {self.forcefield}.')

        os.makedirs(self.side_scaff, exist_ok=True)

        for position in self.mutations:
            chain, resid = position['chain'], position['resid']
            for resname in position['mutations']:
                # Check if self-interaction energy file already exists
                outfile = os.path.join(self.side_scaff, f'{chain}_{resid}_{resname}.csv')
                if os.path.isfile(outfile):
                    logger.info(
                        f'Sidechain-Scaffold interaction energy for residue: {chain}_{resid}_{resname} already computed.')
                    continue
                # Create Molecule object from rotamers of each residue
                # Check if rotamers are computed for all mutations
                try:
                    mol = Molecule(os.path.join(self.rotamer_path, f'{chain}_{resid}', f'{resname}.pdb'))
                except FileNotFoundError:
                    logger.error(f'Missing rotamer for residue: {chain}_{resid}_{resname}.')
                    raise FileNotFoundError(f'Missing rotamer for residue: {chain}_{resid}_{resname}.')

                # get number of rotamers
                nconfs = mol.coords.shape[-1]

                # get paths to mutated scaffolds
                structure_path = os.path.join(self.work_dir, 'scaffold', self.forcefield,
                                              'protein_params', f'{chain}_{resid}_{resname}')

                # Check if directories containing single mutated scaffolds are existing
                if not os.path.isfile(os.path.join(structure_path, 'structure.pdb')):
                    logger.error(f'Missing single mutated structure for mutation: {chain}_{resid}_{resname}.')
                    raise FileNotFoundError(f'Missing single mutated structure for mutation: {chain}_{resid}_{resname}.')

                logger.info(f'Sidechain-Scaffold interaction energy for residue: {chain}_{resid}_{resname} not computed yet.')

                # Read in structure and parameters
                struc, prm = load_ff_parameters(structure_path=structure_path,
                                                forcefield=self.forcefield)
                struc.remove('segid L', _logger=False)

                # Select only the mutated residues without backbone
                variable_selection = f'chain {chain} and resid {resid} and not ({backbone_atoms})'
                # Select all except the flexible side-chains, dont select the backbone because of high Vdw energies
                fixed_selection = f'not segid L and not (chain {chain} and resid {resid} or '
                for res in self.mutations:
                    flex_chain, flex_res = res['chain'], res['resid']
                    if flex_chain != chain or flex_res != resid:
                        # Exclude all other mutated sidechains on different positions, since scored seperately
                        fixed_selection += f'chain {flex_chain} and resid {flex_res} and not ({backbone_atoms}) or '
                fixed_selection = fixed_selection[0:-4] + ')'

                # Generate FFEvaluate object for scoring between sidechain and rest of protein without sidechain
                ffev = FFEvaluate(struc, prm, betweensets=(fixed_selection, variable_selection))

                self_nrgs = np.zeros((nconfs, len(self.energy_terms)))

                with mp.Pool(processes=self.ncpus) as pool:
                    with tqdm(total=nconfs) as pbar:
                        pbar.set_description(f'{chain}_{resid}_{resname}/Scaffold')
                        for index, energy in enumerate(pool.imap(partial(
                                self.calculate_energy,
                                struc=struc,
                                ffev=ffev,
                                res_coords=mol.coords,
                                resid=resid,
                                chain=chain), np.arange(nconfs),
                                chunksize=calculate_chunks(nposes=nconfs, ncpus=self.ncpus))):
                            self_nrgs[index] = energy
                            pbar.update()

                # Save data as csv
                write_energies(outpath=outfile,
                               energies=self_nrgs,
                               energy_terms=self.energy_terms,
                               name_a=resname,
                               nconfs_a=nconfs)

        logger.info('Sidechain-Scaffold calculation was successful.')
