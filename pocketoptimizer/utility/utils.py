import glob
import os
from os.path import dirname
from typing import List, Dict, Tuple, Union, NoReturn

from natsort import natsorted
import itertools
import logging
import parmed
import numpy as np
import pandas as pd
from moleculekit.molecule import Molecule

logger = logging.getLogger(__name__)


class MutationProcessor:

    def __init__(self, structure: str, mutations: List[Dict[str, Union[str, List[str]]]], forcefield: str):
        """
        Constructor Method.

        Class providing functions to process a list of mutations:

            - Remove duplicates
            - Unfold histidine protonation states
            - Verify no mutation is involved in a cystine bridge

        Parameters
        ----------
        structure : str
            Path to the prepared and minimized PDB structure file.
        mutations: list
            List of Dictionaries containing mutations with their corresponding resids and chains
        forcefield: str
            Force field used
        """
        try:
            self.structure = Molecule(structure)
        except FileNotFoundError:
            logger.error(f'Could not find: {structure}.')
            raise FileNotFoundError(f'Could not find: {structure}.')
        self.mutations = mutations
        self.aa = {
            'ALL': ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
                    'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'],
            'AROMATIC': ['PHE', 'TRP', 'TYR'],
            'AMIDE': ['ASN', 'GLN'],
            'ALIPHATIC': ['GLY', 'ALA', 'VAL', 'LEU', 'ILE'],
            'ACIDIC': ['ASP', 'GLU'],
            'BASIC': ['LYS', 'ARG'],
            'HYDRO': ['SER', 'THR'],
            'SULF': ['CYS', 'MET'],
        }
        if forcefield == 'amber_ff14SB':
            self.aa.update({'HIS': ['HID', 'HIE', 'HIP']})
        elif forcefield == 'charmm36':
            self.aa.update({'HIS': ['HSD', 'HSE', 'HSP']})
        self.aa['ALL'].extend(self.aa['HIS'])

    def check_positions(self, check_termini: bool = True) -> NoReturn:
        """
        Checks whether a mutation position exists in the structure,
        is not at the end of a segment and not involved in a disulfide bond
        Removes positions that do not fulfill these criteria

        Parameters
        ----------
        check_termini: bool
            Whether to check for termini positions [default: True]
        """
        updated_mutations = []
        for mutation in self.mutations:
            chain = str(mutation['chain'])
            resid = str(mutation['resid'])
            segid = self.structure.get("segid", sel=f"chain {chain} and resid {resid} and name CA")
            # Check whether there is a C alpha atom belonging to the chain residue identifier combination
            index = self.structure.get('index', f'chain {chain} and resid {resid} and name CA')
            cystine = self.structure.get('resname', sel=f'chain {chain} and resid {resid} and name CA') == 'CYX'

            if not len(index):
                logger.warning(f'Position: {chain}_{resid} is not in the protein.')
                continue
            elif cystine:
                logger.warning(f'Position: {chain}_{resid} was removed because it is involved in a disulfide bond.')
                continue
            elif check_termini and int(resid) == min(self.structure.get('resid', sel=f'segid {segid[0]}')):
                logger.info(f'Position: {chain}_{resid} is the N-terminus of segment: {segid[0]}.')
            elif check_termini and int(resid) == max(self.structure.get('resid', sel=f'segid {segid[0]}')):
                logger.info(f'Position: {chain}_{resid} is the C-terminus of segment: {segid[0]}.')
            else:
                updated_mutations.append(mutation)

        self.mutations = updated_mutations

    def unfold_keywords(self) -> NoReturn:
        """
        - Make mutations upper case
        - Unfold keyword arguments and histidine protonation states
        - Remove non-existing amino acids
        """
        updated_mutations = []

        for i, position in enumerate(self.mutations):
            chain = str(position['chain'])
            resid = str(position['resid'])
            mutations = position['mutations']

            # Make mutations upper case
            mutations = [str(mutation).upper().strip() for mutation in mutations]
            corrected_mutations = []

            # Replace keywords and check if aa is defined
            for mutation in mutations:
                if mutation in self.aa['ALL'] or mutation in self.aa['HIS']:
                    corrected_mutations.append(mutation)
                elif mutation in self.aa.keys():
                    if mutation == 'HIS':
                        logger.warning('Histidine will be unfolded to its three protonation states. '
                                       f'Write {", ".join(self.aa["HIS"])} to define a specific protonation state at N_delta, N_epsilon or N_delta and N_epsilon respectively.')
                    corrected_mutations.extend(self.aa[mutation])
                else:
                    logger.warning(f'Unkown residue: {mutation} has been removed.')

            updated_mutations.append({'mutations': corrected_mutations, 'resid': resid, 'chain': chain})

        self.mutations = updated_mutations

    def merge_duplicates(self) -> NoReturn:
        """
        Merges positions that have been defined multiple times.
        """
        positions = {}
        for i, position in enumerate(self.mutations):
            chain = position['chain']
            resid = position['resid']
            mutations = position['mutations']
            positions.setdefault(f'{chain}_{resid}', mutations).extend(mutations)

        no_duplicates = {}
        for position, mutations in positions.items():
            no_duplicates[position] = list(set(mutations))

        updated_mutations = []
        for position, mutations in no_duplicates.items():
            updated_mutations.append({'mutations': mutations, 'resid': position.split('_')[1], 'chain': position.split('_')[0]})

        self.mutations = updated_mutations

    def check_termini(self) -> List[str]:
        """
        Checks which mutations are intended at the first or last position of a segment

        Returns
        -------

        List containing mutations at termini positions
        """

        termini_positions = []

        for mutation in self.mutations:
            chain = mutation['chain']
            resid = mutation['resid']
            segid = self.structure.get("segid", sel=f"chain {chain} and resid {resid} and name CA")
            if int(resid) == min(self.structure.get('resid', sel=f'segid {segid[0]}')):
                logger.warning(f'Position: {chain}_{resid} is the N-terminus of segment: {segid[0]}.')
                termini_positions.append(f'{chain}_{resid}')
            elif int(resid) == max(self.structure.get('resid', sel=f'segid {segid[0]}')):
                logger.warning(f'Position: {chain}_{resid} is the C-terminus of segment: {segid[0]}.')
                termini_positions.append(f'{chain}_{resid}')

        return termini_positions

    def sort_mutations(self) -> NoReturn:
        """
        First sort mutation lists alphabetically, then sort mutation positions by residue and chain identifier

        """

        mutations_sorted = []
        for position in self.mutations:
            mutations_sorted.append({'mutations': natsorted(position['mutations']), 'resid': position['resid'], 'chain': position['chain']})

        # Sorting self.mutations first by chain then by resid
        mutations_sorted = natsorted(mutations_sorted, key=lambda k: (k['chain'], k['resid']))
        self.mutations = mutations_sorted

    def process_mutations(self, check_termini: bool = True) -> List[Dict[str, Union[str, List[str]]]]:
        """
        Process mutations

        Parameters
        ----------
        check_termini: bool
            Whether to check for termini positions [default: True]

        Return
        ------
        Process list of mutations
        """
        # Check that position in protein and not involved in cystine bridge
        self.check_positions(check_termini=check_termini)
        # Replace keyword arguments and check for amino acids
        self.unfold_keywords()
        # merge duplicates
        self.merge_duplicates()
        # Sort mutations alphabetically
        self.sort_mutations()

        return self.mutations


def create_pairs(mutations: List[Dict[str, Union[str, List[str]]]]) -> List[List[str]]:
    """
    Creates pairs of mutations based on a 2D-array.

    Parameters
    ----------
    mutations: list
        List of Dictionaries containing mutations with their corresponding resids and chains

    Returns
    -------
    List of Lists containing naturally sorted pairwise mutations
    [['A_44_THR', 'B_54_ARG'], ['A_44_THR', 'B_54_ILE'], ['A_44_THR', 'B_57_TRP'],
    ['B_54_TRP', 'B_57_LEU'], ['B_54_TRP', 'B_57_VAL']]
    """
    position = []
    for pos in mutations:
        # All mutation positions
        position.append(['_'.join(x) for x in itertools.product([pos['chain'] + "_" + pos['resid']], pos['mutations'])])
    pair_permutations = []
    # Create all possible pairwise permutations of two residues at different positions
    for i, j in itertools.combinations(position, 2):
        pair_permutations.extend(x for x in itertools.product(i, j))

    # Natural sort permutation tuples to prevent getting wrong directory

    sorted_pair_permutations = []
    for pair in pair_permutations:
        sorted_pair_permutations.append(natsorted(pair))

    return sorted_pair_permutations


def fix_parameters(parameterfile: str) -> str:
    """
    Somehow the CHARMM .prm files have all the stuff
    commented. This function uncomments it.

    Parameters
    ----------
    parameterfile: str
        Path to a directory containing a structure.prm file.

    Returns
    -------
    Path to the fixed file.
    """
    new_prm_file = os.path.join(dirname(parameterfile), 'structure.prm')
    with open(parameterfile, 'r') as f:
        lines = f.readlines()

    with open(new_prm_file, 'w') as f:
        for line in lines:
            line = line.replace('!MASS', 'MASS')
            line = line.replace('!ATOMS', 'ATOMS')
            f.write(line)
    return new_prm_file


def load_ff_parameters(structure_path: str, forcefield: str) -> Tuple[Molecule, Union[parmed.charmm.CharmmParameterSet, parmed.amber.AmberParameterSet]]:
    """
    Convenience function to load Amber or CHARMM parameter files.

    Parameters
    ----------
    structure_path: str
        Path to a directory containing structure.prmtop file for amber or structure.prm and structure.psf files for charmm, as well as structure.pdb file
    forcefield: str
        Forcefield for which the system was built

    Returns
    -------
    Parameters and :class: moleculekit.molecule.Molecule object of the structure build for the respective force field calculations
    """
    import warnings
    warnings.filterwarnings("ignore")

    if forcefield.startswith('charmm'):
        prmfiles = glob.glob(os.path.join(structure_path, 'topologies', '*.rtf'))
        prmfiles.append(os.path.join(structure_path, 'structure.prm'))
        prm = parmed.charmm.CharmmParameterSet(*prmfiles)
        struc = Molecule(os.path.join(structure_path, 'structure.psf'))
        struc.read(os.path.join(structure_path, 'structure.pdb'))
    elif forcefield.startswith('amber'):
        _struct = parmed.amber.AmberParm(os.path.join(structure_path, 'structure.prmtop'), )
        prm = parmed.amber.AmberParameterSet.from_structure(_struct)
        struc = Molecule(os.path.join(structure_path, 'structure.prmtop'))
        struc.read(os.path.join(structure_path, 'structure.pdb'))
    else:
        logger.error(f'Force field: {forcefield} not supported.')
        raise NotImplementedError(f'Force field: {forcefield} not supported.')
    return struc, prm


def write_energies(outpath: str, energies: np.ndarray, energy_terms: List[str], name_a: str, nconfs_a: int, name_b: str = '', nconfs_b: int = 0) -> NoReturn:
    """
    Writes energy values stored in a numpy array into a .csv file

    Parameters
    ----------
    outpath: str
        Path to write csv file to
    energies: np.ndarray
        Energy array
    energy_terms: list
        List of energy terms to write as columns
    name_a: str
        Name of residue/ligand
    nconfs_a: int
        Number of conformations/poses
    name_b: str
        Name of residue b [default: '']
    nconfs_b: int
        Number of conformations [default: 0]

    Returns
    -------
    """
    # Round energy values to four decimal places
    energies = np.around(energies, 4)

    # Convert numpy array to dataframe
    df = pd.DataFrame(energies)

    # Create row names
    rows = [f'{name_a}_{conf}' for conf in np.arange(nconfs_a)]
    df.insert(loc=0, column=name_a, value=rows)

    # Create column names
    columns = [name_a]
    if not nconfs_b:
        for energy_term in energy_terms:
            columns.append(energy_term)
    else:
        for conf in np.arange(nconfs_b):
            for energy_term in energy_terms:
                columns.append(f'{name_b}_{conf}_{energy_term}')
    df.columns = columns

    df.to_csv(outpath, header=True, index=False, sep='\t')


def calculate_chunks(nposes: int, ncpus: int) -> int:
    """
    Calculates the size of the chunks for multiprocessing by dividing the number of
    poses by the number of available CPUs times 4

    Parameters
    ----------
    nposes: int
        Number of poses
    ncpus: int
        Number of available CPUs

    Returns
    -------
    Size of each chunk
    """

    chunksize, extra = divmod(nposes, ncpus * 4)
    if extra:
        chunksize += 1
    return chunksize


class DotDict(dict):
    """
    Example:
    m = Map({'first_name': 'Eduardo'}, last_name='Pool', age=24, sports=['Soccer'])
    """

    def __init__(self, *args, **kwargs):
        super(DotDict, self).__init__(*args, **kwargs)
        for arg in args:
            if isinstance(arg, dict):
                for k, v in arg.items():
                    self[k] = v

        if kwargs:
            for k, v in kwargs.items():
                self[k] = v

    def __getattr__(self, attr):
        return self.get(attr)

    def __setattr__(self, key, value):
        self.__setitem__(key, value)

    def __setitem__(self, key, value):
        super(DotDict, self).__setitem__(key, value)
        self.__dict__.update({key: value})

    def __delattr__(self, item):
        self.__delitem__(item)

    def __delitem__(self, key):
        super(DotDict, self).__delitem__(key)
        del self.__dict__[key]
