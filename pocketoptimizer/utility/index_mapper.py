# In order to type hint a class object in a class method in python > 3.7
from __future__ import annotations
import os
from typing import List, Dict, Tuple, Union, Generator, NoReturn
import logging
from natsort import natsorted
from moleculekit.molecule import Molecule

from pocketoptimizer.ui import DesignPipeline
from pocketoptimizer.utility.utils import Storer

logger = logging.getLogger(__name__)


class IndexMapper(Storer):
    """
    Class used to record for a design the names of the design positions
    (backbone positions, ligand, waters, ...), the names of
    residues allowed at these positions, and the number of conformers
    of each residues. It maps the positions, residues and conformers
    to indices in raw energy files and solver
    input files, so they can be mapped to each other and retrieved.
    """

    def __init__(self, **kwargs):
        """
        Constructor.
        Creates empty member lists and dictionaries.
        Can be filled by from_index_file or from_conformer_files

        self._conformer_counts: Contains the number of conformers for each position
                                {'B_41':{'ALA':1, 'ASP':13, ...}, 'B_44':{..},
                                'ligand':{'ligand':506}, .. }

        self._scaling_factors: Scaling Factor, {'ligand': 10.0, 'metal':10.0, ..}
                               If no factor is included for a given design position, 1.0 is assumed

        self._dummy_energy: This energy value (in kcal/mol) is set as the pairwise energy of
                            conformer pairs that are excluded
        """
        super().__init__(**kwargs)
        self._conformer_counts = {}
        self._scaling_factors = {}
        self._dummy_energy = 1e20

    def _read_index_file(self, filename: str) -> NoReturn:
        """
        Reads an index file.

        Parameters
        ----------
        filename: str
            Path of the index file
        """
        self._conformer_counts = {}
        self._scaling_factors = {}
        section = ""
        with open(filename) as f:
            for line in f:
                if line.startswith('['):
                    section = line.strip()
                elif line.isspace():
                    continue
                elif section == '[POSITIONS]':
                    pos, dummy = line.split()
                    self._conformer_counts[pos] = {}
                elif section == '[CONFORMERS]':
                    p, dummy, res, dummy = line.split()
                    if res in self._conformer_counts[p]:
                        self._conformer_counts[p][res] += 1
                    else:
                        self._conformer_counts[p][res] = 1
                elif section == '[SCALING]':
                    pos, sf = line.split()
                    self._scaling_factors[pos] = float(sf)

    @classmethod
    def from_index_file(cls, filename: str, design_pipeline: DesignPipeline) -> IndexMapper:
        """
        Initialize from an index file.

        Parameters
        ----------
        filename: str
            Path of the index file
        design_pipeline: :class: DesignPipeline
            DesignPipeline object for initialization
        """
        im = cls(design_pipeline=design_pipeline)
        im._read_index_file(filename)
        return im

    @classmethod
    def from_conformer_files(cls, design_pipeline: DesignPipeline) -> IndexMapper:
        """
        Initializes IndexMapper object from side chain rotamer and ligand conformer files.

        Parameters
        ----------
        design_pipeline: :class: DesignPipeline
            DesignPipeline object for initialization

        Returns
        --------
        :class:IndexMapper object
        """
        im = cls(design_pipeline=design_pipeline)
        im._count_conformers()
        return im

    def _count_conformers(self) -> NoReturn:
        """
        Initializes IndexMapper object with _conformers_count attribute. Dictionary containing
        Positions as keys ('A_44') and dictionaries as values. These dictionaries contain resnames for positions
        as keys and corresponding conformer counts as values.
        """
        self._conformer_counts = {}
        self._scaling_factors = {}
        for position in self.mutations:
            chain, resid = position['chain'], position['resid']
            self._conformer_counts[f'{chain}_{resid}'] = {}
            for resname in position['mutations']:
                # Try to open rotamer file and count the number of Models contained
                try:
                    with open(os.path.join(self.rotamer_path, f'{chain}_{resid}', f'{resname}.pdb'), 'r') as rotamer_file:
                        for line in rotamer_file:
                            if line.startswith('MODEL'):
                                self._conformer_counts[f'{chain}_{resid}'][resname] = \
                                    self._conformer_counts[f'{chain}_{resid}'].get(resname, 0) + 1
                except:
                    logger.error(f'Missing rotamer for residue: {chain}_{resid}_{resname}.')
                    raise FileNotFoundError(f'Missing rotamer for residue: {chain}_{resid}_{resname}.')
        try:
            ligand_poses = Molecule(self.ligand_poses_pdb)
            ligand_poses.read(self.ligand_poses_xtc)
        except:
            logger.error(f'{self.ligand_poses_pdb} not found.')
            raise FileNotFoundError(f'{self.ligand_poses_pdb} not found.')

        # Count the number of frames, which corresponds to the number of ligand poses
        self._conformer_counts['ligand'] = {}
        self._conformer_counts['ligand']['ligand'] = ligand_poses.coords.shape[-1]

    def write_index_file(self, filename: str) -> NoReturn:
        """
        Writes an index file that maps the indices of the ligand poses and
        the ids of the flexible backbone positions and conformers
        to the indices used for the kingsford file. Needed to process a
        kingsford solution output.

        Parameters
        ----------
        filename: str
            Path of the index file to be created
        """

        # Write the position indices
        with open(filename, 'w') as f:
            f.write('[POSITIONS]\n')
            for index, position in enumerate(natsorted(self._conformer_counts)):
                f.write(f'{position} {index}\n')
            f.write('\n')
            f.write('[CONFORMERS]\n')
            for position in natsorted(self._conformer_counts):
                position_confs = 0
                for res in natsorted(self._conformer_counts[position]):
                    for i in range(self._conformer_counts[position][res]):
                        f.write(f'{position} {position_confs} {res} {i}\n')
                        position_confs += 1
            f.write('\n')
            f.write('[SCALING]\n')
            for pose in natsorted(self._scaling_factors):
                f.write(f'{pose} {self._scaling_factors[pose]}\n')
            f.write('\n')

    def set_exclude_res(self, exclude: List[Dict[str, Union[str, List[str]]]]) -> NoReturn:
        """
        Excludes residues from the design procedure.

        Parameters
        ----------
        exclude: list
            Mutations at certain positions to be exlcuded from the design
            [{'mutations': ['ALA', 'ASN', 'GLU'], 'resid': '23', 'chain': 'A'},
            {'mutations': ['ALA', 'SER'], 'resid': '27', 'chain': 'A'}]
        """
        for residue in exclude:
            _resname = f'{residue["chain"]}_{residue["resid"]}'
            if _resname in self._conformer_counts:
                if 'ALL' in residue['mutations']:
                    self._conformer_counts.pop(_resname, None)
                else:
                    for aa in residue['mutations']:
                        if aa in self._conformer_counts[_resname]:
                            self._conformer_counts[_resname].pop(aa, None)
                    if self._conformer_counts[_resname] == {}:
                        self._conformer_counts.pop(_resname, None)

    def set_dummy_energy(self, energy: float) -> NoReturn:
        """
        Set the dummy energy value.

        Parameters
        ----------
        energy: float
            Energy value
        """
        self._dummy_energy = str(energy)

    def get_pos_count(self) -> int:
        """
        Counts the number of design positions.

        Returns
        -------
        Number of design positions (i.e. backbone positions, ligand(s), water(s), metal(s))
        """
        return len(self._conformer_counts)

    def get_conf_count_for_pos(self, position: str) -> int:
        """
        Gets number of conformers for a position.

        Parameters
        ----------
        position: str
            Design position (e.g. 'B_41', 'ligand', 'water1')

        Returns
        -------
        The total number of conformers (for all residue types) at given position.
        """

        c = 0
        for dummy_res, i in self._conformer_counts[position].items():
            c += i
        return c

    def get_conf_count_for_res(self, pos: str, res: str) -> int:
        """
        Gets number of conformers for a residue at a position.

        Parameters
        ----------
        pos: str
            Design position (e.g. 'B_41', 'ligand', 'water1', 'metal', ...)
        res: str
            Residue type in three letter code ('ARG')

        Returns
        -------
        The number of conformers for one residue type at given position.
        """
        return self._conformer_counts[pos][res]

    def get_resnames(self, position: str) -> List[str]:
        """
        Gets list of resnames for a position.

        Parameters
        ----------
        position: str
            Design position (e.g. 'B_41', 'ligand', 'water1')

        Returns
        -------
        The names of the residues for which there are conformers at a given position
        """
        return natsorted(self._conformer_counts[position])

    def get_pos_names(self) -> List[str]:
        """
        Gets naturally sorted list of all design positions.

        Returns
        -------
        The names of the design positions ordered naturally
        """
        return natsorted(self._conformer_counts)

    def get_restype_and_index(self, pos: str, index: int) -> Tuple[str, int]:
        """
        Get the residue type and index of a conformer at a position.

        Parameters
        ----------
        pos: str
            Design position (e.g. 'B_41', 'ligand', 'water1', ...)
        index: int
            Global index of a conformer at a position (e.g. taken from the kingsford energy file)

        Returns
        -------
        A tuple (res, index) where res is the residue name and index the index of the conformer of that residue
        """
        sorty = natsorted(self._conformer_counts[pos])
        base = 0
        for ir, r in enumerate(sorty):
            if ir == len(sorty) - 1:
                return r, index - base
            elif base + self._conformer_counts[pos][r] > index >= base:
                return r, index - base
            else:
                base += self._conformer_counts[pos][r]
        logger.error(f'No residue and/or index found for {pos}, {index}.')
        raise 'No residue and/or index found for ' + pos + ', ' + index

    def get_restype(self, pos: str, index: int) -> str:
        """
        Get residue name for a global index and a position.

        Parameters
        ----------
        pos: str
            Design position (e.g. 'B_41', 'ligand', 'water1', ...)
        index: int
            Global index of a conformer at a position (e.g. taken from the kingsford energy file)

        Returns
        -------
        The residue name
        """
        base = 0
        for r in natsorted(self._conformer_counts[pos]):
            if index > self._conformer_counts[pos][r] + base:
                base += self._conformer_counts[pos][r]
            else:
                return r

    def get_pos_name(self, index: int) -> str:
        """
        Get name of a design position from index.

        Parameters
        ----------
        index: int
            Position index from the kingsford input file

        Returns
        -------
        The name of that design position, (e.g. 'B_41', 'ligand', 'water1', ...)
        """
        return natsorted(self._conformer_counts)[index]

    def get_pos_index(self, pos: str) -> int:
        """
        Get a position index from the name of a design position.

        Parameters
        ----------
        pos: int
            Position name, (e.g. 'B_41', 'ligand', 'water1', ...)

        Returns
        -------
        Index of the position
        """
        return natsorted(self._conformer_counts).index(pos)

    def get_pairs(self) -> Generator[Tuple[str], None, None]:
        """
        Get the names of two design positions following each other.

        Returns
        -------
        Tuple containing the two design positions
        """
        for i, pose_1 in enumerate(self.get_pos_names()):
            for pose_2 in self.get_pos_names()[i + 1:]:
                yield (pose_1, pose_2)

    def set_scaling_factor(self, pos: str, factor: float) -> NoReturn:
        """
        Set scaling factor for energies involving one position.

        Parameters
        ----------
        pos: str
            Design position (e.g. 'B_41', 'ligand', 'water1', ...)
        factor: float
            New scaling factor
        """
        self._scaling_factors[pos] = factor

    def set_scaling_factors(self, pos_to_factor: Dict[str, float]) -> NoReturn:
        """
        Set scaling factor for energies involving several position.

        Parameters
        ----------
        pos_to_factor: dict
            Position names as keys and scaling factors as values
        """
        for pose, factor in pos_to_factor.items():
            self._scaling_factors[pose] = factor

    def get_scaling_factor(self, pos: str) -> float:
        """
        Get the energy scaling factor for a position

        Parameters
        ----------
        pos: str
            Design position (e.g. 'B_41', 'ligand', 'water1', ...)
        """
        if pos in self._scaling_factors:
            return self._scaling_factors[pos]
        else:
            return 1.0
