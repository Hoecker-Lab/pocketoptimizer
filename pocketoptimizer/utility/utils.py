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
from rdkit import Chem

logging.root.handlers = []
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - [%(levelname)s] - %(message)s",
    handlers=[
        logging.FileHandler(os.environ.get('POCKETOPTIMIZER_LOGFILE')),
        logging.StreamHandler()
    ]
)

logger = logging.getLogger('pocketoptimizer.utilities')


class MutationProcessor:

    def __init__(self, scaffold: str, mutations: List[Dict[str, Union[str, List[str]]]]):
        """
        Constructor Method.

        Class providing functions to process a list of mutations:

            - Remove duplicates
            - Unfold histidine protonation states
            - Verify no mutation is involved in a cystine bridge

        Parameters
        ----------
        scaffold : str
            Path to the prepared and minimized PDB structure file.
        mutations: list
            List of Dictionaries containing mutations with their corresponding resids and chains
        """
        try:
            self.structure = Molecule(scaffold)
        except FileNotFoundError:
            logger.error(f'Could not find: {scaffold}.')
            raise FileNotFoundError(f'Could not find: {scaffold}.')
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
            'HIS': ['HID', 'HIE', 'HIP']
        }
        self.mutations = mutations

    def check_positions(self) -> NoReturn:
        """
        Checks whether a mutation position exists in the structure,
        is not at the end of a segment and not involved in a disulfide bond
        Removes positions that do not fulfill these criteria
        """
        updated_mutations = []
        for mutation in self.mutations:
            chain = str(mutation['chain'])
            resid = str(mutation['resid'])
            segid = self.structure.get("segid", sel=f"chain {chain} and resid {resid} and name CA")
            # Check whether there is a C alpha atom belonging to the chain residue identifier combination
            index = self.structure.get('index', f'chain {chain} and resid {resid} and name CA')
            cystine = self.structure.get('resname', sel=f'chain {chain} and resid {resid} and name CA') == 'CYX'

            include = True

            if not len(index):
                include = False
                logger.warning(f'Position: {chain}_{resid} is not in the protein.')
            if cystine:
                logger.warning(f'Position: {chain}_{resid} was removed because it is involved in a disulfide bond.')
                include = False
            if int(resid) == min(self.structure.get('resid', sel=f'segid {segid[0]}')):
                logger.warning(f'Position: {chain}_{resid} is the N-terminus of segment: {segid[0]}.')
                include = False
            elif int(resid) == max(self.structure.get('resid', sel=f'segid {segid[0]}')):
                logger.warning(f'Position: {chain}_{resid} is the C-terminus of segment: {segid[0]}.')
                include = False

            if include:
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
            mutations = [str(mutation).upper() for mutation in mutations]
            corrected_mutations = []

            # Replace keywords and check if aa is defined
            for mutation in mutations:
                if mutation in self.aa['ALL'] or mutation in self.aa['HIS']:
                    corrected_mutations.append(mutation)
                elif mutation in self.aa.keys():
                    if mutation == 'HIS':
                        logger.warning('Histidine will be unfolded to its three protonation states. '
                                       'Write HID, HIE, HIP to define a specific protonation state at N_delta, N_epsilon or N_delta and N_epsilon respectively.')
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

    def check_termini(self) -> Dict[str, List[str]]:
        """
        Checks which mutations are intended at the first or last position of a segment

        Returns
        -------

        Dictionary containing lists of mutations sorted after N-and C-terminus
        """

        terminus_positions = {}

        for mutation in self.mutations:
            chain = mutation['chain']
            resid = mutation['resid']
            segid = self.structure.get("segid", sel=f"chain {chain} and resid {resid} and name CA")
            if int(resid) == min(self.structure.get('resid', sel=f'segid {segid[0]}')):
                logger.warning(f'Position: {chain}_{resid} is the N-terminus of segment: {segid[0]}.')
                terminus_positions.setdefault('N-terminus', [f'{chain}_{resid}']).append(f'{chain}_{resid}')
            elif int(resid) == max(self.structure.get('resid', sel=f'segid {segid[0]}')):
                logger.warning(f'Position: {chain}_{resid} is the C-terminus of segment: {segid[0]}.')
                terminus_positions.setdefault('C-terminus', [f'{chain}_{resid}']).append(f'{chain}_{resid}')

        return terminus_positions

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

    def process_mutations(self) -> List[Dict[str, Union[str, List[str]]]]:
        """
        Process mutations
        """
        # Check that position in protein and not involved in cystine bridge
        self.check_positions()
        # Replace keyword arguments and check for amino acids
        self.unfold_keywords()
        # merge duplicates
        self.merge_duplicates()
        # Sort mutations alphabetically
        self.sort_mutations()

        return self.mutations

# EDIT 03/2021: If TRP is used as an input ligand in .mol file type, try converting it to .mol2 with obabel
# and remove atom names to make it not be recognized as a protein.
def is_trp(ligand: str) -> bool:
    """
    Check if an input ligand is TRP.

    Parameters
    ----------
    ligand: str
        Path to ligand

    Returns
    -------
    True if ligand is TRP, otherwise False
    """
    with open(ligand) as rf:
        from collections import defaultdict
        atoms = defaultdict(int)
        found_first_atm = False
        for line in rf:
            if len(line.split()) > 7:
                found_first_atm = True
                atoms[line.split()[3]] += 1
            elif found_first_atm:
                break
        if "C" in atoms and "N" in atoms and "O" in atoms:
            if atoms["C"] == 11 and atoms["N"] == 2 and atoms["O"] == 2:
                return True
    return False


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
    # Round energy values to three decimal places
    energies = np.around(energies, 3)

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
