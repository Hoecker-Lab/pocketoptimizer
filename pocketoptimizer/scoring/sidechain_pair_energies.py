import itertools
import multiprocessing as mp
import os
from functools import partial
from typing import List, Tuple, Dict, Union, NoReturn
import logging
import numpy as np
from ffevaluation.ffevaluate import FFEvaluate
from moleculekit.molecule import Molecule
from tqdm.auto import tqdm

logger = logging.getLogger(__name__)


def calculate_energy(ids: Tuple[int], structure: Molecule, ffev: FFEvaluate, res_a_coords: np.ndarray, res_b_coords: np.ndarray,\
                     resid_a: str, resid_b: str, chain_a: str, chain_b: str) -> Tuple[np.float]:
    """
    Calculating the energy between two rotamers at two different positions

    Parameters
    ----------
    ids: Tuple
        Indices of current rotamers
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
    structure.set('coords', res_a_coords[:, :, ids[0]], f'chain {chain_a} and resid {resid_a}')
    structure.set('coords', res_b_coords[:, :, ids[1]], f'chain {chain_b} and resid {resid_b}')
    energies = ffev.calculateEnergies(structure.coords)
    return energies['vdw'], energies['elec'] * 0.01


def calculate_pairs(work_dir: str, rotamer_path: str, mutations: List[Dict[str, Union[str, List[str]]]],
                    forcefield: str, ncpus: int = 1) -> NoReturn:
    """
    Wrapper function to calculate sidechain sidechain energies using molecular mechanics force field energy computation

    Parameters
    ----------
    work_dir: str
        Path to working directory
    rotamer_path: str
        Path containing sidechain rotamer .pdb files
    mutations: List[Dict]
        List of dictionaries defining the mutations at different positions
    forcefield: str
        Forcefield to use for scoring
    ncpus: int
        Number of CPUs to use for sidechain scaffold scoring [default: 1]
    """
    from pocketoptimizer.utility.utils import load_ff_parameters, create_pairs, write_energies, calculate_chunks
    from pocketoptimizer.utility.molecule_types import backbone_atoms

    logger.info(f'Compute energies using forcefield: {forcefield}.')

    output_dir = os.path.join(work_dir, 'energies', f'{forcefield}_{rotamer_path.split("/")[-1]}', f'sidechain_sidechain_{forcefield}')
    os.makedirs(output_dir, exist_ok=True)

    sorted_pair_permutations = create_pairs(mutations)

    for mutation_a, mutation_b in sorted_pair_permutations:
        # Get name of each residue from each pair
        chain_a, resid_a = mutation_a.split('_')[0:2]
        chain_b, resid_b = mutation_b.split('_')[0:2]
        chain_a_resid_a = '_'.join(mutation_a.split('_')[0:2])
        chain_b_resid_b = '_'.join(mutation_b.split('_')[0:2])
        resname_a = mutation_a.split('_')[2]
        resname_b = mutation_b.split('_')[2]
        # Check if pairwise-interaction energy files already exists
        outfile = os.path.join(output_dir,f'{chain_a}_{resid_a}_{resname_a}_{chain_b}_{resid_b}_{resname_b}.csv')
        if os.path.isfile(outfile):
            logger.info(f'Sidechain-Sidechain-Interaction-Energies for Residue-Pair: {chain_a}_{resid_a}_{resname_a}/{chain_b}_{resid_b}_{resname_b} already computed.')
            continue
        # Create Molecule object from rotamers of each residue
        # Check if rotamers are computed for all mutations
        try:
            res_a = Molecule(os.path.join(rotamer_path, chain_a_resid_a, f'{resname_a}.pdb'))
        except FileNotFoundError:
            logger.error(f'Missing rotamer for residue: {chain_a_resid_a}_{resname_a}.')
            raise FileNotFoundError(f'Missing rotamer for residue: {chain_a_resid_a}_{resname_a}.')
        try:
            res_b = Molecule(os.path.join(rotamer_path, chain_b_resid_b, f'{resname_b}.pdb'))
        except FileNotFoundError:
            logger.error(f'Missing rotamer for residue: {chain_b_resid_b}_{resname_b}.')
            raise FileNotFoundError(f'Missing rotamer for residue: {chain_b_resid_b}_{resname_b}.')

        # Create output path from residues of pair
        structure_path = os.path.join(work_dir, 'scaffold', forcefield, 'protein_params',
                                      f'{chain_a}_{resid_a}_{resname_a}_'
                                      f'{chain_b}_{resid_b}_{resname_b}')
        # Check if directories containing pairwise mutated scaffolds are existing
        if not os.path.isfile(os.path.join(structure_path, 'structure.pdb')):
            logger.error(f'Missing pairwise mutated structure for mutation combination: {chain_a}_{resid_a}_{resname_a}_{chain_b}_{resid_b}_{resname_b}.')
            raise FileNotFoundError(f'Missing pairwise mutated structure for mutation combination: {chain_a}_{resid_a}_{resname_a}_{chain_b}_{resid_b}_{resname_b}.')

        logger.info(
            f'Sidechain-Sidechain-Interaction-Energies for Residue-Pair: {chain_a}_{resid_a}_{resname_a}/{chain_b}_{resid_b}_{resname_b} not computed yet.')

        struc, prm = load_ff_parameters(structure_path=structure_path, forcefield=forcefield)

        sel_res_a = f'chain {chain_a} and resid {resid_a} and not ({backbone_atoms})'
        sel_res_b = f'chain {chain_b} and resid {resid_b} and not ({backbone_atoms})'

        # Create FFEvaluate object between residue a and b without backbone
        ffev = FFEvaluate(struc, prm, betweensets=(sel_res_a, sel_res_b))

        # List of all possible rotamer index combinations for two mutations at two different positions [(0,0), (0,1)..]
        nconfs_a = res_a.coords.shape[-1]
        nconfs_b = res_b.coords.shape[-1]

        batch = list(itertools.product(*[np.arange(nconfs_a), np.arange(nconfs_b)]))

        # Create array (nconfs * nconfs * 2 (vdw and es))
        pair_nrgs = np.zeros((nconfs_a, nconfs_b * 2))

        with mp.Pool(processes=ncpus) as pool:
            with tqdm(total=len(batch)) as pbar:
                pbar.set_description(
                    f'{chain_a}_{resid_a}_{resname_a}/'
                    f'{chain_b}_{resid_b}_{resname_b}')
                for index, energy in enumerate(pool.imap(partial(
                        calculate_energy,
                        structure=struc.copy(),
                        ffev=ffev,
                        res_a_coords=res_a.coords,
                        res_b_coords=res_b.coords,
                        resid_a=resid_a,
                        resid_b=resid_b,
                        chain_a=chain_a,
                        chain_b=chain_b),
                        batch,
                        chunksize=calculate_chunks(nposes=len(batch), ncpus=ncpus))):
                    pair_nrgs[batch[index][0], batch[index][1]*2:batch[index][1]*2 + 2] = energy
                    pbar.update()

                # Save data as csv
                write_energies(outpath=outfile,
                               energies=pair_nrgs,
                               energy_terms=['VdW', 'ES'],
                               name_a=resname_a,
                               nconfs_a=nconfs_a,
                               name_b=resname_b,
                               nconfs_b=nconfs_b)
    logger.info('Sidechain-Pair calculation was successful.')
