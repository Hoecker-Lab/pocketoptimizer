import itertools
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


class SidechainPairScorer(Storer):
    """
    Scorer class for sidechain/sidechain interactions
    Utilizing molecular dynamics binding energy related computations
    Requires built complexes
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.energy_terms = ['Vdw', 'Elec']

    def calculate_energy(self, ids: Tuple[int], struc: Molecule, ffev: FFEvaluate,
                         res_a_coords: np.ndarray, res_b_coords: np.ndarray,
                         resid_a: str, resid_b: str, chain_a: str, chain_b: str) -> np.ndarray:
        """
        Calculating the energy between two rotamers at two different positions

        Parameters
        ----------
        ids: tuple
            Tuple of rotamer ids
        structure: :class: moleculekit.molecule.Molecule
            Object containing all rotamers as frames
        ffev: :class: ffevaluation.ffevaluate.FFEvaluate
            Object created for respective force field
            from parameter file between sets fixed_selection and variable_selection,
            where energies can be calculated from
        res_a_coords: np.ndarray
            Coordinates of the rotamers of residue a
        res_b_coords: np.ndarray
            Coordinates of the rotamers of residue b
        resid_a: str
            Residue ID of the mutated residue a
        resid_b: str
            Residue ID of the mutated residue b
        chain_a: str
            Chain ID of the mutated residue a
        chain_b: str
            Chain ID of the mutated residue b

        Returns
        -------
        Tuple with energy components

        """
        struc.set('coords', res_a_coords[:, :, ids[0]], f'chain {chain_a} and resid {resid_a}')
        struc.set('coords', res_b_coords[:, :, ids[1]], f'chain {chain_b} and resid {resid_b}')
        energies = ffev.calculateEnergies(struc.coords)
        return np.array([energies['vdw'], energies['elec'] * 0.01])

    def calculate_pairs(self) -> NoReturn:
        """
        Wrapper function to calculate sidechain sidechain energies using molecular mechanics force field energy computation
        """

        from pocketoptimizer.utility.utils import load_ff_parameters, create_pairs, write_energies, calculate_chunks
        from pocketoptimizer.utility.molecule_types import backbone_atoms

        logger.info(f'Compute energies using forcefield: {self.forcefield}.')

        os.makedirs(self.side_side, exist_ok=True)

        sorted_pair_permutations = create_pairs(self.mutations)

        for mutation_a, mutation_b in sorted_pair_permutations:
            # Get name of each residue from each pair
            chain_a, resid_a = mutation_a.split('_')[0:2]
            chain_b, resid_b = mutation_b.split('_')[0:2]
            chain_a_resid_a = '_'.join(mutation_a.split('_')[0:2])
            chain_b_resid_b = '_'.join(mutation_b.split('_')[0:2])
            resname_a = mutation_a.split('_')[2]
            resname_b = mutation_b.split('_')[2]
            # Check if pairwise-interaction energy files already exists
            outfile = os.path.join(self.side_side, f'{chain_a}_{resid_a}_{resname_a}_{chain_b}_{resid_b}_{resname_b}.csv')
            if os.path.isfile(outfile):
                logger.info(f'Sidechain-Sidechain interaction energy for residue pair: '
                            f'{chain_a}_{resid_a}_{resname_a}/{chain_b}_{resid_b}_{resname_b} already computed.')
                continue
            # Create Molecule object from rotamers of each residue
            # Check if rotamers are computed for all mutations
            try:
                res_a = Molecule(os.path.join(self.rotamer_path, chain_a_resid_a, f'{resname_a}.pdb'))
                res_a.remove(backbone_atoms, _logger=False)
            except FileNotFoundError:
                logger.error(f'Missing rotamer for residue: {chain_a_resid_a}_{resname_a}.')
                raise FileNotFoundError(f'Missing rotamer for residue: {chain_a_resid_a}_{resname_a}.')
            try:
                res_b = Molecule(os.path.join(self.rotamer_path, chain_b_resid_b, f'{resname_b}.pdb'))
                res_b.remove(backbone_atoms, _logger=False)
            except FileNotFoundError:
                logger.error(f'Missing rotamer for residue: {chain_b_resid_b}_{resname_b}.')
                raise FileNotFoundError(f'Missing rotamer for residue: {chain_b_resid_b}_{resname_b}.')

            # Create output path from residues of pair
            structure_path = os.path.join(self.work_dir, 'scaffold', self.forcefield, 'protein_params',
                                          f'{chain_a}_{resid_a}_{resname_a}_'
                                          f'{chain_b}_{resid_b}_{resname_b}')
            # Check if directories containing pairwise mutated scaffolds are existing
            if not os.path.isfile(os.path.join(structure_path, 'structure.pdb')):
                logger.error(f'Missing pairwise mutated structure for mutation combination: '
                             f'{chain_a}_{resid_a}_{resname_a}_{chain_b}_{resid_b}_{resname_b}.')
                raise FileNotFoundError(f'Missing pairwise mutated structure for mutation combination: '
                                        f'{chain_a}_{resid_a}_{resname_a}_{chain_b}_{resid_b}_{resname_b}.')

            logger.info(
                f'Sidechain-Sidechain interaction energy for residue pair: {chain_a}_{resid_a}_{resname_a}/{chain_b}_{resid_b}_{resname_b} not computed yet.')

            struc, prm = load_ff_parameters(structure_path=structure_path, forcefield=self.forcefield)

            sel_res_a = f'chain {chain_a} and resid {resid_a} and not ({backbone_atoms})'
            sel_res_b = f'chain {chain_b} and resid {resid_b} and not ({backbone_atoms})'

            sidechains = struc.copy()
            sidechains.filter(sel=f'((chain {chain_a} and resid {resid_a}) or (chain {chain_b} and resid {resid_b})) and not ({backbone_atoms})', _logger=False)

            # Create FFEvaluate object between residue a and b without backbone
            ffev = FFEvaluate(sidechains, prm, betweensets=(sel_res_a, sel_res_b))

            # List of all possible rotamer index combinations for two mutations at two different positions [(0,0), (0,1)..]
            nconfs_a = res_a.coords.shape[-1]
            nconfs_b = res_b.coords.shape[-1]

            batch = list(itertools.product(*[np.arange(nconfs_a), np.arange(nconfs_b)]))

            # Create array (nconfs * nconfs * 2 (vdw and es))
            pair_nrgs = np.zeros((nconfs_a, nconfs_b * len(self.energy_terms)))

            with mp.Pool(processes=self.ncpus) as pool:
                with tqdm(total=len(batch)) as pbar:
                    pbar.set_description(
                        f'{chain_a}_{resid_a}_{resname_a}/'
                        f'{chain_b}_{resid_b}_{resname_b}')
                    for index, energy in enumerate(pool.imap(partial(
                            self.calculate_energy,
                            struc=sidechains,
                            ffev=ffev,
                            res_a_coords=res_a.coords,
                            res_b_coords=res_b.coords,
                            resid_a=resid_a,
                            resid_b=resid_b,
                            chain_a=chain_a,
                            chain_b=chain_b),
                            batch,
                            chunksize=calculate_chunks(nposes=len(batch), ncpus=self.ncpus))):
                        pair_nrgs[batch[index][0], batch[index][1]*2:batch[index][1]*2 + 2] = energy
                        pbar.update()

                    # Save data as csv
                    write_energies(outpath=outfile,
                                   energies=pair_nrgs,
                                   energy_terms=self.energy_terms,
                                   name_a=resname_a,
                                   nconfs_a=nconfs_a,
                                   name_b=resname_b,
                                   nconfs_b=nconfs_b)

        logger.info('Sidechain-Pair calculation was successful.')
