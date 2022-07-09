import multiprocessing as mp
import os
from functools import partial
from typing import Tuple, List, Dict, Union, NoReturn
import logging
import numpy as np
from ffevaluation.ffevaluate import FFEvaluate
from moleculekit.molecule import Molecule
from tqdm.auto import tqdm

logging.root.handlers = []
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - [%(levelname)s] - %(message)s",
    handlers=[
        logging.FileHandler(os.environ.get('POCKETOPTIMIZER_LOGFILE')),
        logging.StreamHandler()
    ]
)

logger = logging.getLogger('pocketoptimizer.scoring.sidechain_scaffold_energies')


def calculate_energy(id: int, structure: Molecule, ffev: FFEvaluate, res_coords: np.ndarray, resid: str, chain: str) -> Tuple[np.float]:
    """
    Calculating the energy between a rotamer and the fixed protein scaffold

    Parameters
    ----------
    id: int
        Index of current rotamer
    structure: :class: moleculekit.molecule.Molecule
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
    Tuple with energy components
    """

    # Set coordinates to current rotamer
    structure.set('coords', res_coords[:, :, id], f'chain {chain} and resid {resid}')
    energies = ffev.calculateEnergies(structure.coords)
    return energies['vdw'], energies['elec']


def calculate_scaffold(work_dir: str, rotamer_path: str, mutations: List[Dict[str, Union[str, List[str]]]], forcefield: str, ncpus: int = 1) -> None or NoReturn:
    """
    Wrapper function to calculate sidechain scaffold energies using molecular mechanics force field energy computation

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
    from pocketoptimizer.utility.utils import load_ff_parameters, write_energies, calculate_chunks

    logger.info(f'Compute energies using forcefield: {forcefield}.')

    output_dir = os.path.join(work_dir, 'energies', f'{forcefield}_{rotamer_path.split("/")[-1]}', f'sidechain_scaffold_{forcefield}')
    os.makedirs(output_dir, exist_ok=True)

    for position in mutations:
        chain, resid = position['chain'], position['resid']
        for resname in position['mutations']:
            # Check if self-interaction energy file already exists
            outfile = os.path.join(output_dir, f'{chain}_{resid}_{resname}.csv')
            if os.path.isfile(outfile):
                logger.info(
                    f'Sidechain-Scaffold-Interaction-Energies for Residue: {chain}_{resid}_{resname} already computed.')
                continue
            # Create Molecule object from rotamers of each residue
            # Check if rotamers are computed for all mutations
            try:
                mol = Molecule(os.path.join(rotamer_path, f'{chain}_{resid}', f'{resname}.pdb'))
            except FileNotFoundError:
                logger.error(f'Missing rotamer for residue: {chain}_{resid}_{resname}.')
                raise FileNotFoundError(f'Missing rotamer for residue: {chain}_{resid}_{resname}.')

            # get number of rotamers
            nconfs = mol.coords.shape[-1]

            # get paths to mutated scaffolds
            structure_path = os.path.join(work_dir, 'scaffold', forcefield, 'protein_params', f'{chain}_{resid}_{resname}')

            # Check if directories containing single mutated scaffolds are existing
            if not os.path.isfile(os.path.join(structure_path, 'structure.pdb')):
                logger.error(f'Missing single mutated structure for mutation: {chain}_{resid}_{resname}.')
                raise FileNotFoundError(f'Missing single mutated structure for mutation: {chain}_{resid}_{resname}.')

            logger.info(f'Sidechain-Scaffold-Interaction-Energies for Residue: {chain}_{resid}_{resname} not computed yet.')

            # Read in structure and parameters
            struc, prm = load_ff_parameters(structure_path=structure_path, forcefield=forcefield)

            # Select only the mutated residues without backbone
            variable_selection = f'chain {chain} and resid {resid} and sidechain'
            # Select all except the mutated residues, dont select the backbone because of high Vdw energies
            fixed_selection = f'protein and not (chain {chain} and resid {resid} or '
            for res in mutations:
                flex_chain, flex_res = res['chain'], res['resid']
                if flex_chain != chain or flex_res != resid:
                    # Exclude all other mutated sidechains on different positions, since scored seperately
                        fixed_selection += f'chain {flex_chain} and resid {flex_res} and sidechain or '
            fixed_selection = fixed_selection[0:-4] + ')'

            # Generate FFEvaluate object for scoring between sidechain and rest of protein without sidechain
            ffev = FFEvaluate(struc, prm, betweensets=(fixed_selection, variable_selection))

            self_nrgs = np.zeros((nconfs, 2))

            with mp.Pool(processes=ncpus) as pool:
                with tqdm(total=nconfs) as pbar:
                    pbar.set_description(f'{chain}_{resid}_{resname}/Scaffold')
                    for index, energy in enumerate(pool.imap(partial(
                            calculate_energy,
                            structure=struc.copy(),
                            ffev=ffev,
                            res_coords=mol.coords,
                            resid=resid,
                            chain=chain), np.arange(nconfs),
                            chunksize=calculate_chunks(nposes=nconfs, ncpus=ncpus))):
                        self_nrgs[index] = energy
                        pbar.update()

            # Save data as csv
            write_energies(outpath=outfile,
                           energies=self_nrgs,
                           energy_terms=['VdW', 'ES'],
                           name_a=resname,
                           nconfs_a=nconfs)

    logger.info('Sidechain-Scaffold calculation was successful.')
