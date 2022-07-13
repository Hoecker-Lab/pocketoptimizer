import copy

from tqdm.auto import tqdm
import os
import glob
from shutil import which, rmtree
from typing import List, Dict, Union, NoReturn
import logging

logger = logging.getLogger('pocketoptimizer.ui')

# API
class DesignPipeline:
    """This class acts as a UI for accessing key features of PocketOptimizer.

    It contains all necessary steps required for a full design run.
    Everything is build around a user defined work directory in which a
    specific directory tree will be created during the design run.
    The steps should be executed in the following order:

    1. Structure Preparation:
        - parameterize ligand
        - prepare_protein
        - prepare_mutants
    2. Modelling of flexibility
        - prepare_lig_conformers
        - sample_sidechain_rotamers
        - sample_lig_poses
    3. Energy Calculation
        - calculate_energies
    4. Calculate the best solutions
    """

    def __init__(self, work_dir: str, forcefield: str, ph: float = 7.2, peptide: bool = False, ncpus: int = 1):
        """
        Constructor method

        Parameters
        ----------
        work_dir: str
            Working directory, all files will be created in this directory
        forcefield: str
            Force field used for sampling procedures and energy calculation steps [options: 'amber_ff14SB', 'charmm36']
        pH: float
            pH value that will be used to adjust the protonation states [default: 7.2]
        peptide: bool
            Whether the ligand is a peptide [default: False]
        ncpus:
            Number of CPUs to use [default: 1]
        """

        self.work_dir = work_dir
        os.chdir(self.work_dir)
        # Set environment variable to standard path
        os.environ["POCKETOPTIMIZER_LOGFILE"] = os.path.join(self.work_dir, 'pocketoptimizer.log')
        from pocketoptimizer.path import path

        # Set environment variables for MATCH
        os.environ["MATCH"] = os.path.join(path(), '..', 'MATCH_RELEASE', 'MATCH')
        os.environ["PerlChemistry"] = os.path.join(path(), '..', 'MATCH_RELEASE', 'PerlChemistry')

        logging.root.handlers = []
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s - %(name)s - [%(levelname)s] - %(message)s",
            handlers=[
                logging.FileHandler(os.environ.get('POCKETOPTIMIZER_LOGFILE')),
                logging.StreamHandler()
            ]
        )

        def create_settings_py(work_dir: str) -> NoReturn:
            """
            Creates a settings.py file in the working directory which contains paths to all binary files,
            the paths are accessible as attributes of the Settings class.

            Parameters
            ----------

            work_dir: Path to working directory
            """
            import tempfile

            if not os.path.isdir(work_dir):
                logger.error('Working directory does not exist.')
                raise FileNotFoundError('Working directory does not exist.')

            tmp_dir = tempfile.gettempdir()
            pocketoptimizer_logfile = os.path.join(work_dir, 'pocketoptimizer.log')

            match_path = os.path.join(path(), '..', 'MATCH_RELEASE', 'MATCH', 'scripts', 'MATCH.pl')
            solver_path = path(bin_file='sontag_solver')
            psfgen_path = path(bin_file='psfgen')

            obabel_path = which('obabel')
            if not obabel_path:
                obabel_path = '#'
                logger.warning("Obabel not found make sure it is installed and in your path.")
            antechamber_path = which('antechamber')
            if not antechamber_path:
                antechamber_path = '#'
                logger.warning("Antechamber not found make sure ambertools/ambermini is installed and in your path.")
            tleap_path = which('tleap')
            if not tleap_path:
                tleap_path = '#'
                logger.warning("Tleap not found make sure ambertools/ambermini is installed and in your path.")
            parmchk2_path = which('parmchk2')
            if not parmchk2_path:
                parmchk2_path = '#'
                logger.warning("Parmchk2 not found make sure it is installed and in your path.")
            smina_path = which('smina')
            if not smina_path:
                smina_path = '#'
                logger.warning("Smina not found make sure it is installed and in your path.")

            if os.path.isfile(os.path.join(work_dir, 'settings', 'settings.py')):
                logger.warning('Settings.py already exists, procede with caution.')

            else:
                settings_folder = os.path.join(work_dir, 'settings')
                os.makedirs(settings_folder, exist_ok=True)
                os.mknod(os.path.join(settings_folder, '__init__.py'))
                with open(os.path.join(settings_folder, 'settings.py'), 'w') as settings_file:
                    settings_file.write(f"OBABEL_BIN = '{obabel_path}'\n")
                    settings_file.write(f"MATCH_PERL = '{match_path}'\n")
                    settings_file.write(f"ANTECHAMBER_BIN = '{antechamber_path}'\n")
                    settings_file.write(f"PARMCHK2_BIN = '{parmchk2_path}'\n")
                    settings_file.write(f"TLEAP_BIN = '{tleap_path}'\n")
                    settings_file.write(f"PSFGEN_BIN = '{psfgen_path}'\n")
                    settings_file.write(f"SMINA_BIN = '{smina_path}'\n")
                    settings_file.write(f"SOLVER_BIN = '{solver_path}'\n")
                    settings_file.write(f"POCKETOPTIMIZER_LOGFILE = '{pocketoptimizer_logfile}'\n")
                    settings_file.write(f"TMP_DIR = '{tmp_dir}'\n")
                    settings_file.write('\nclass Settings:'
                                        '\n\n    def __init__(self):'
                                        '\n        self.OBABEL_BIN = OBABEL_BIN'
                                        '\n        self.MATCH_PERL = MATCH_PERL'
                                        '\n        self.ANTECHAMBER_BIN = ANTECHAMBER_BIN'
                                        '\n        self.PARMCHK2_BIN = PARMCHK2_BIN'
                                        '\n        self.TLEAP_BIN = TLEAP_BIN'
                                        '\n        self.PSFGEN_BIN = PSFGEN_BIN'
                                        '\n        self.SMINA_BIN = SMINA_BIN'
                                        '\n        self.SOLVER_BIN = SOLVER_BIN'
                                        '\n        self.POCKETOPTIMIZER_LOGFILE = POCKETOPTIMIZER_LOGFILE'
                                        '\n        self.TMP_DIR = TMP_DIR')

        create_settings_py(self.work_dir)
        import sys
        import multiprocessing as mp
        if not self.work_dir in sys.path:
            sys.path.insert(0, self.work_dir)
        self._update_settings()
        logger.info(f'Logging to: {self.settings.POCKETOPTIMIZER_LOGFILE}')

        self.built_scaffold = ''
        self.prepared_protein = ''
        self.peptide = peptide
        self.ligand_protonated = ''
        self.built_ligand_params = None
        self.ligand_conformers = ''
        self.ligand_poses = ''
        self.mutations = []
        self.scoring = {
            'smina': ['vina', 'vinardo', 'dkoes_scoring', 'ad4_scoring'],
            'ff': ['amber_ff14SB', 'charmm36'],
        }
        self.scorer = ''
        self.forcefield = ''
        self._set_ff(forcefield)
        self.rotamer_path = ''
        self.ph = ph
        self.temperature = 300.0
        self.ncpus = ncpus
        if self.ncpus > mp.cpu_count():
            logger.warning(f'More CPUs defined than available, setting maximum of {mp.cpu_count()} CPUs.')
            self.ncpus = mp.cpu_count()
        os.environ['NUMEXPR_MAX_THREADS'] = str(self.ncpus)

    def _update_settings(self):
        """
        Reimports settings file to update binary paths, sets environment variable to logfile path and updates logger.

        Returns:
        --------

        Updated :class: settings.Settings object containing the paths currently set in settings.py
        """

        if '#' in open(os.path.join(self.work_dir, 'settings', 'settings.py'), 'r').read():
            logger.warning('Missing binary paths in settings file, adjust manually.')

        from importlib import reload
        import settings.settings as settings
        settings = reload(settings)
        self.settings = settings.Settings()

        # Set environment variable to path defined in settingsfile
        os.environ["POCKETOPTIMIZER_LOGFILE"] = self.settings.POCKETOPTIMIZER_LOGFILE
        logging.root.handlers = []
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s - %(name)s - [%(levelname)s] - %(message)s",
            handlers=[
                logging.FileHandler(os.environ.get('POCKETOPTIMIZER_LOGFILE')),
                logging.StreamHandler()
            ]
        )
        if not os.path.isdir(self.settings.TMP_DIR):
            logger.info(f'Creating temporary directory: {self.settings.TMP_DIR}')
            os.makedirs(self.settings.TMP_DIR, exist_ok=True)

    def _set_ff(self, forcefield: str) -> NoReturn:
        """
        Sets protein and ligand force fields used for all sampling procedures and energy calculations.

        Parameters
        ----------
        forcefield: str
            Force field. Options are 'amber_ff14SB' and 'charmm36' (Default: 'amber_ff14SB').
         """
        from pocketoptimizer.utility.utils import DotDict
        if forcefield in self.scoring['ff']:
            self.forcefield = forcefield
            self.built_scaffold = os.path.join(self.work_dir, 'scaffold', self.forcefield, 'scaffold.pdb')
            self.prepared_protein = os.path.join(self.work_dir, 'scaffold', self.forcefield, 'protein_preparation', 'prepared_scaffold.pdb')
            if not self.peptide:
                self.ligand_protonated = os.path.join(self.work_dir, 'ligand', self.forcefield, 'ligand.mol2')
            else:
                self.ligand_protonated = os.path.join(self.work_dir, 'ligand', self.forcefield, 'ligand.pdb')
            self.built_ligand_params = DotDict({'params_folder': os.path.join(self.work_dir, 'ligand', self.forcefield, 'params'),
                                                'mol2': [os.path.join(self.work_dir, 'ligand', self.forcefield, 'params','ligand.mol2')],
                                                'frcmod': [os.path.join(self.work_dir, 'ligand', self.forcefield, 'params','ligand.frcmod')],
                                                'prm': [os.path.join(self.work_dir, 'ligand', self.forcefield, 'params','ligand.prm')],
                                                'rtf': [os.path.join(self.work_dir, 'ligand', self.forcefield, 'params','ligand.rtf')]})
            self.ligand_conformers = os.path.join(self.work_dir, 'ligand', self.forcefield, 'conformers', 'ligand_confs.pdb')
            self.ligand_poses = os.path.join(self.work_dir, 'ligand', self.forcefield, 'poses', 'ligand_poses.pdb')
        else:
            logger.error(f'{forcefield} not supported. Use {self.scoring["ff"][0]} or {self.scoring["ff"][1]}.')
            raise ValueError(f'{forcefield} not supported. Use {self.scoring["ff"][0]} or {self.scoring["ff"][1]}.')

    def _set_rot_lib(self, rot_lib: str = 'dunbrack') -> NoReturn:
        """
        Sets the rotamer path, based on the rotamer libary

        Parameters
        ----------
        rot_lib: str
            Rotamer library used for rotamer sampling.
            Options are either: 'dunbrack' or 'cmlib' [default: 'dunbrack']
        """
        if rot_lib not in ['dunbrack', 'cmlib']:
            logger.error('Rotamer library not implemented.')
            raise NotImplementedError('Rotamer library not implemented.')
        self.rotamer_path = os.path.join(self.work_dir, 'scaffold', self.forcefield,
                                             'rotamers', rot_lib)

    def set_temperature(self, temperature: float) -> NoReturn:
        """
        Sets the temperature used in the design

        Parameters
        ----------
        temperature: temperature value to set

        Returns
        -------

        """
        try:
            self.temperature = float(temperature)
        except ValueError:
            logger.error('Temperature must be a number.')
            raise ValueError('Temperature must be a number.')

    def set_mutations(self, mutations: List[Dict[str, Union[str, List[str]]]]) -> NoReturn:
        """
        Sets the mutations that will be used to prepare the structures for side chain rotamer sampling and scoring.
        Does also support the following keywords: 'ALL', 'AROMATIC', 'AMIDE', 'ALIPHATIC', 'ACIDIC', 'BASIC',
        'HYDRO', 'SULF'.

        'ALL': ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE',
                'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'],
        'AROMATIC': ['PHE', 'TRP', 'TYR', 'HIS'],
        'AMIDE': ['ASN', 'GLN'],
        'ALIPHATIC': ['GLY', 'ALA', 'VAL', 'LEU', 'ILE'],
        'ACIDIC': ['ASP', 'GLU'],
        'BASIC': ['LYS', 'ARG', 'HIS'],
        'HYDRO': ['SER', 'THR'],
        'SULF': ['CYS', 'MET']

        Parameters
        ----------
        mutations: list
            List of Dictionaries containing mutations with their corresponding resids and chains
        """
        from pocketoptimizer.utility.utils import MutationProcessor

        logger.info('If design positions are removed or added a new design run should be started.')

        mutation_processor = MutationProcessor(structure=self.prepared_protein,
                                               mutations=mutations)
        self.mutations = mutation_processor.process_mutations()

    def parameterize_ligand(self, input_ligand: str) -> NoReturn:
        """
        Interface for Ligand parameterization.
        Adds hydrogen atoms to the input molecule, produces GAFF2 or CGenFF_2b6 parameters.

        Parameters
        ----------
        input_ligand: str
            Path to a ligand file
        """
        from pocketoptimizer.preparation.structure_building import SystemBuilder

        self._update_settings()

        if not os.path.isfile(self.ligand_protonated):

            system = SystemBuilder(work_dir=self.work_dir,
                                   structure=os.path.join(self.work_dir, input_ligand),
                                   forcefield=self.forcefield,
                                   ligand_params=self.built_ligand_params,
                                   tleap=self.settings.TLEAP_BIN,
                                   psfgen=self.settings.PSFGEN_BIN,
                                   obabel=self.settings.OBABEL_BIN,
                                   antechamber=self.settings.ANTECHAMBER_BIN,
                                   parmchk2=self.settings.PARMCHK2_BIN,
                                   match=self.settings.MATCH_PERL)

            system.parameterize_ligand(ph=self.ph)

        else:
            logger.info('Ligand is already parametrized.')

    def prepare_lig_conformers(self, method: str = 'genetic', nconfs: int = 100,
                               score: str = 'rmsd', rcutoff: float = 0.5, ecutoff: float = 50.0) -> NoReturn:
        """
        Samples ligand conformers to account for ligand flexibility during the design process.
        Writes conformers into ligand_confs.pdb into ligand/conformers in project directory.

        Parameters
        ----------
        method: str
            Method that will be used to sample conformers. Implemented options are 'confab' and 'genetic' [default: genetic].
        nconfs: int
            Maximum number of conformers to sample (default 100). Final number can be lower
            due to pruning or diversity filtering.
        score: str
            Scoring option to decide how to filter the generated conformers [default: rmsd].
        rcutoff: float
            RMSD cutoff (default 0.5 Angstrom). Only used together with the 'confab' method. [default: 0.5 Angstrom]
        ecutoff: float
            Energy cutoff (default 50.0 kcal/mol). Only used together with the 'confab' method. [default: 50 kcal/mol]
        """
        from pocketoptimizer.sampling import conformer_generator_obabel

        self._update_settings()
        os.makedirs(os.path.join(self.work_dir, 'ligand', self.forcefield, 'conformers'), exist_ok=True)

        if method == 'genetic' or method == 'confab':
            conformer_generator_obabel.conformer_generator(
                obabel_path=self.settings.OBABEL_BIN, infile=self.ligand_protonated, conf_file_name=self.ligand_conformers,
                method=method, nconfs=nconfs, score=score, rcutoff=rcutoff, ecutoff=ecutoff)
        else:
            logger.error('Conformer generation method not defined! Try genetic or confab.')
            raise ValueError('Conformer generation method not defined! Try genetic or confab.')

    def prepare_peptide_conformers(self, positions: List[str], library: str = 'dunbrack',
                                   nrg_thresh: float = 100.0, expand: List[str] = ['chi1', 'chi2'],
                                   dunbrack_filter_thresh: float = 0.01):
        """
        positions: list
            List of residue IDs to sample rotamers for
        library: str
            Library to use for selecting rotamers, either setting coordinates from pdb after superimposing (cmlib)
            or setting dihedral angles from dunbrack [default: 'dunbrack']
        nrg_thresh: float
            Threshold value used to filter peptide conformations [default: 100.0 kcal/mol]
        expand: list
            List of which chi angles to expand, [default: ['chi1', 'chi2']]
        dunbrack_filter_thresh: float
            Only rotamers with an occurence above this value are selected from the library [default: 0.01 (1%)]
        """
        from pocketoptimizer.sampling.peptide_sampler import PeptideSampler

        self._update_settings()

        if nrg_thresh < 0:
            logger.error('Threshold needs to be positive.')
            return

        _positions = []
        for resid in set(positions):
            _positions.append({'chain': 'L', 'resid': resid})

        os.makedirs(os.path.join(self.work_dir, 'ligand', self.forcefield, 'conformers'), exist_ok=True)

        peptide_sampler = PeptideSampler(work_dir=self.work_dir,
                                         positions=_positions,
                                         forcefield=self.forcefield,
                                         params_folder=self.built_ligand_params.params_folder,
                                         library=library,
                                         )

        peptide_sampler.conformer_sampling(nrg_thresh=nrg_thresh,
                                           dunbrack_filter_thresh=dunbrack_filter_thresh,
                                           expand=expand,
                                           ncpus=self.ncpus)
        os.chdir(self.work_dir)

    def prepare_protein(self, protein_structure: str, keep_chains: List = [], backbone_restraint: bool = True,
                        discard_mols: List[Dict[str, str]] = [], peptide_structure: str = None,
                        peptide_mutations: List[Dict[str, str]] = [], cuda: bool = False) -> NoReturn or None:
            """
            Protonates and cleans the protein structure followed by a subsequent minimization step.

            First step in the design pipeline is the preparation of the scaffold.
            This includes adding hydrogens to the structure and protonating polar
            sidechains accordingly.
            Furthermore, the structure gets cleaned of unwanted chains, waters, ions etc.
            If the structures already exist, the procedure will be skipped.
            After protein preparation the structure will be minimized with the ligand
            inside the binding pocket.

            Parameters
            ----------
            protein_structure: str
                Path to the protein PDB file.
            keep_chains: list
                Protein chain which will be extracted for the design (['A', 'B'])
            backbone_restraint: bool
                Restraints the backbone during minimization. [default: True]
            discard_mols: list
                List of dictionaries containing chain and resid of the molecules to discard. In the case of amino acid ligands,
                these have to be manually discarded with this option. [default: None]
            peptide: str
                Path to peptide that should be prepared
            peptide_mutations: list
                allows to set peptide mutations in the following way [{'resid': '6', 'mutation': 'ALA'}] [default: None]
            cuda: bool
                If the minimization should be performed on a GPU. [default: False]
            """
            from pocketoptimizer.preparation import minimize_structure
            from pocketoptimizer.preparation.structure_building import SystemBuilder

            self._update_settings()

            logger.info('Start Protein Preparation.')

            if self.peptide:
                if not peptide_structure:
                    logger.error('No input peptide structure defined.')
                    raise RuntimeError('No input peptide structure defined.')

                else:
                    _peptide_mutations = []
                    for mutation in peptide_mutations:
                        _peptide_mutations.append({'chain': 'L',
                                                   'resid': mutation['resid'],
                                                   'mutation': mutation['mutation']})

                    system = SystemBuilder(work_dir=self.work_dir,
                                           structure=os.path.join(self.work_dir, peptide_structure),
                                           forcefield=self.forcefield,
                                           ligand_params=self.built_ligand_params,
                                           mutations=_peptide_mutations,
                                           peptide=self.peptide,
                                           tleap=self.settings.TLEAP_BIN,
                                           psfgen=self.settings.PSFGEN_BIN)
                    system.prepare_peptide(prepared_peptide=self.ligand_protonated, pH=self.ph)
                    # Build force field parameters for peptide
                    system.build_peptide()

            if not type(keep_chains) == list:
                logger.error('Define chains in a list.')
                raise TypeError('Define chains in a list.')

            if not type(discard_mols) == list:
                logger.error('Define molecules in a list of dictionaries.')
                raise TypeError('Define molecules in a list of dictionaries.')

            if not type(peptide_mutations) == list:
                logger.error('Define mutations in a list of dictionaries.')
                raise TypeError('Define mutations in a list of dictionaries.')

            system = SystemBuilder(work_dir=self.work_dir,
                                   structure=os.path.join(self.work_dir, protein_structure),
                                   forcefield=self.forcefield,
                                   peptide=self.peptide,
                                   ligand_params=self.built_ligand_params,
                                   tleap=self.settings.TLEAP_BIN,
                                   psfgen=self.settings.PSFGEN_BIN)

            os.makedirs(os.path.join(self.work_dir, 'scaffold', self.forcefield, 'protein_preparation'), exist_ok=True)
            system.prepare_protein(prepared_protein=self.prepared_protein,
                                   keep_chains=keep_chains,
                                   pH=self.ph,
                                   discard_mols=discard_mols)

            # Check if there is a built scaffold
            if not os.path.isfile(self.built_scaffold):
                # Minimize the structure with ligand in pocket
                logger.info('Building complex before minimization.')
                # Build the complex before minimization
                system.build_complex(ligand=self.ligand_protonated)
                input_path = os.path.join(self.work_dir, 'scaffold', self.forcefield, 'protein_params', 'native_complex')

                if not self.peptide:
                    minimize_structure.minimize_structure(
                        structure_path=input_path,
                        forcefield=self.forcefield,
                        output_pdb=self.built_scaffold,
                        cuda=cuda,
                        restraint_bb=backbone_restraint,
                        temperature=self.temperature)
                else:
                    minimize_structure.minimize_structure(
                        structure_path=input_path,
                        forcefield=self.forcefield,
                        output_pdb=self.built_scaffold,
                        cuda=cuda,
                        restraint_bb=backbone_restraint,
                        minimize_ligand=True,
                        output_ligand=self.ligand_protonated,
                        temperature=self.temperature)
                logger.info('Your protein was successfully minimized and can be used for design now.')

            else:
                logger.info('Your scaffold is already minimized.')

            if os.path.exists(os.path.join(self.work_dir, 'scaffold', 'conversion_table.txt')):
                logger.info(f'Please change your residue ID numbers according to {os.path.join(self.work_dir, "scaffold", "conversion_table.txt")} and start from set_mutations() again.')

            logger.info('Protein preparation finished.')

    def prepare_mutants(self, sampling_pocket: str = 'GLY') -> NoReturn:
        """
        Prepares protein mutants for computations with a specific force field.

        Builds all user defined protein mutation variants of single and pairwise mutated scaffolds
        in order to compute self an pairwise interaction energies for packing and binding contributions.
        Also builds a full glycin pockets to sample ligand poses and rotamers.

        Writes structures, parameters etc. into the '/my_project/scaffold/protein_params' project directory.

        sampling_pocket: str
            What type of sampling pockets to build for rotamer and ligand pose sampling procedures
            In theory every pocket is possible, but only glycine and alanine are meanigful [default: 'GLY']
        """
        from pocketoptimizer.preparation.structure_building import SystemBuilder
        from pocketoptimizer.preparation.aacodes import aa

        self._update_settings()

        if not self.mutations:
            logger.error('No mutations have been defined. Please define at least one mutation.')
            raise AttributeError('No mutations have been defined. Please define at least one mutation.')

        logger.info('Start building mutated protein scaffold variants.')

        if sampling_pocket not in aa:
            logger.error('Amino acid unknown.')
            raise ValueError('Amino acid unknown.')
        logger.info(f'Build {sampling_pocket} sampling pockets.')

        system = SystemBuilder(work_dir=self.work_dir,
                              structure=self.built_scaffold,
                              mutations=self.mutations,
                              forcefield=self.forcefield,
                              ligand_params=self.built_ligand_params,
                              peptide=self.peptide,
                              tleap=self.settings.TLEAP_BIN,
                              psfgen=self.settings.PSFGEN_BIN)
        # Build all single and pairwise mutated scaffolds and all sampling pockets for sampling rotamers and ligand poses
        system.build_complex(ligand=self.ligand_protonated, sampling_pocket=sampling_pocket)
        logger.info('Scaffold building done.')

    def sample_sidechain_rotamers(self, library: str = 'dunbrack', vdw_filter_thresh: float = 100.0,
                                  dunbrack_filter_thresh: float = 0.01, expand: List[str] = ['chi1', 'chi2']) -> NoReturn:
        """Uses the dunbrack or cm_lib rotamer library to compute possible rotamers of defined residues.
        Energies for rotamer pruning are calculated using ffevaluate with the amber_ff14SB or charmm36 force field

        Writes residue wise rotamer.pdb's into scaffold/rotamers in the form of:
        'scaffold/rotamers/A_171/ASP.pdb'

        Parameters
        ----------
        library: str
            Rotamer library, options are: 'cm_lib' or 'dunbrack', 'dunbrack' is recommended [default: 'dunbrack']
        vdw_filter_thresh: float
            Threshold value used to filter clashing rotamers [default: 100.0 kcal/mol]
        dunbrack_filter_thresh: float
            Filter threshold, rotamers having probability of occurence lower than filter threshold will
            be pruned if their rotameric mod does occur more than once, value should be between 0 and 1
            (-1: no pruning, 1: pruning of all rotamers with duplicate rotamer modes) [default: 0.01]
        expand: list
            List of which chi angles to expand, [default: ['chi1', 'chi2']]
        """

        from pocketoptimizer.sampling import sidechain_rotamers_ffev

        self._update_settings()

        if not self.mutations:
            logger.error('No mutations have been defined. Please define at least one mutation.')
            raise AttributeError('No mutations have been defined. Please define at least one mutation.')
        # set rotamer path for energy calculations
        self._set_rot_lib(library)

        if vdw_filter_thresh < 0:
            logger.error('Threshold needs to be positive.')
            return

        rotamer_sampler = sidechain_rotamers_ffev.FFRotamerSampler(
            work_dir=self.work_dir,
            mutations=copy.deepcopy(self.mutations),
            forcefield=self.forcefield,
            rot_path=self.rotamer_path,
            library=library,
            tmp=self.settings.TMP_DIR)

        rotamer_sampler.rotamer_sampling(
            dunbrack_filter_thresh=dunbrack_filter_thresh,
            expand=expand,
            vdw_filter_thresh=vdw_filter_thresh,
            ncpus=self.ncpus)
        os.chdir(self.work_dir)

    def sample_lig_poses(self, method: str = 'grid', grid: Dict[str, List[float]] = None,
                         vdw_filter_thresh: float = 100.0, max_poses: int = 10000) -> NoReturn:
        """Calculate possible ligand poses inside a glycin/alanin pocket.

        Filters all poses above a certain vdW energy threshold to remove
        poses that clash with the backbone.
        Writes ligand pose into ligand/ligand_poses.pdb and additional poses
        as a trajectory file in ligand/ligand_poses.xtc.

        Parameters
        ----------
        method: str
            Method to sample the ligand. Options are 'grid' and 'random' [default: 'grid'].
        grid: dict
            Dictionary containing grid parameters for sampling poses.
            Has the following structure:
             {'trans':[DISTANCE, STEPSIZE],
              'rot':[ANGLE, STEPSIZE]}
        vdw_filter_thresh: int
            Energy threshold in kcal/mol [default: 100 kcal/mol]
        max_poses: int
            Maximum amount of poses to produce [default: 10000].
        """
        from pocketoptimizer.sampling import ligand_poses

        if method == 'grid' and not grid:
            logger.error('Define grid, if grid method should be used.')
            raise ValueError('Define grid, if grid method should be used.')

        if vdw_filter_thresh < 0:
            logger.error('Threshold needs to be positive.')
            return

        sampling_scaffold_ligand = os.path.join(self.work_dir, 'scaffold', self.forcefield, 'protein_params',
                                                'ligand_sampling_pocket', 'structure.pdb')

        if os.path.isfile(sampling_scaffold_ligand):
            if not os.path.isfile(self.ligand_poses):
                logger.info('Sample ligand poses.')
                sampler = ligand_poses.PoseSampler(work_dir=self.work_dir,
                                                   ligand=self.ligand_conformers,
                                                   forcefield=self.forcefield,
                                                   max_poses=max_poses,
                                                   sampling_method=method,
                                                   peptide=self.peptide)

                sampler.create_ligand_ensemble(grid=grid, vdw_filter_thresh=vdw_filter_thresh, ncpus=self.ncpus)
            else:
                logger.info('Ligand poses are already sampled.')
        else:
            logger.error('Ligand pose sampling scaffold has not been prepared yet.')

    def calculate_energies(self, scoring: str = 'vina') -> NoReturn:
        """Computes self and pairwise energies of all ligand and scaffold components.

        The energy calculation procedure consists of three parts.
        1. Sidechain vs. Scaffold energies
        2. Sidechain pair energies
        3. Ligand vs. scaffold/sidechain energies

        The former two are related to force field computations while the latter is done using a
        specific scoring function.
        The output energies/scores will be written into .csv files accessible in the
        /work_dir/energies/ directory.

        Parameters
        ----------
        scoring: str
            Method that will be used for scoring [default: 'vina'].
            The following methods are available: vina, vinardo, dkoes_scoring, ad4_scoring, charmm36, amber_ff14SB, nnscore, rfscore, plecscore, cnnscore
        """
        from pocketoptimizer.scoring import sidechain_scaffold_energies, sidechain_pair_energies, smina_scorer, ff_scorer
        os.chdir(self.work_dir)

        self._update_settings()

        if not self.mutations:
            logger.error('No mutations have been defined. Please define at least one mutation.')
            raise AttributeError('No mutations have been defined. Please define at least one mutation.')
        if not self.rotamer_path:
            logger.error('No rotamer library has been defined. Sample side chain rotamers first.')
            raise AttributeError('No rotamer library has been defined. Sample side chain rotamers first.')

        if scoring in self.scoring['smina'] or scoring in self.scoring['ff']:
            self.scorer = scoring
        else:
            logger.error('Scorer not implemented.')
            raise NotImplementedError('Scorer not implemented.')

        logger.info('Start energy calculation.')
        logger.info(f'Using {self.ncpus} CPUs for multiprocessing.')

        logger.info('Calculate Sidechain-Scaffold Energies.')
        sidechain_scaffold_energies.calculate_scaffold(
            work_dir=self.work_dir,
            rotamer_path=self.rotamer_path,
            mutations=copy.deepcopy(self.mutations),
            forcefield=self.forcefield,
            ncpus=self.ncpus

        )
        logger.info('Calculate Sidechain-Pair Energies.')
        sidechain_pair_energies.calculate_pairs(
            work_dir=self.work_dir,
            rotamer_path=self.rotamer_path,
            mutations=copy.deepcopy(self.mutations),
            forcefield=self.forcefield,
            ncpus=self.ncpus

            )
        logger.info('Calculate Ligand-Scaffold/Sidechain-Interaction-Energies.')
        if self.scorer in self.scoring['smina']:
            vs = smina_scorer.SminaScorer(
                project_path=self.work_dir,
                scaffold=self.built_scaffold,
                ligand=self.ligand_protonated,
                ligand_poses=self.ligand_poses,
                scoring_method=self.scorer,
                forcefield=self.forcefield,
                rotamer_path=self.rotamer_path,
                smina_path=self.settings.SMINA_BIN,
                peptide=self.peptide,
                tmp_dir=self.settings.TMP_DIR,
                obabel=self.settings.OBABEL_BIN,
            )

            vs.run_smina_scorer(
                mutations=copy.deepcopy(self.mutations),
                ncpus=self.ncpus

            )
        elif self.scorer in self.scoring['ff']:

            ffs = ff_scorer.FFScorer(
                project_path=self.work_dir,
                ligand=self.ligand_protonated,
                ligand_poses=self.ligand_poses,
                forcefield=self.scorer,
                rotamer_path=self.rotamer_path
            )
            ffs.run_ff_scorer(
                mutations=copy.deepcopy(self.mutations),
                ncpus=self.ncpus
            )

        else:
            logger.error('Scorer not implemented.')
            raise NotImplementedError('Scorer not implemented.')

        os.chdir(self.work_dir)
        logger.info('Energy calculation was successful.')

    def design(self, num_solutions: int, ligand_scaling: int = 1) -> NoReturn:
        """
        Computes defined number of solutions that minimise the total energy and creates output files.
        In three step process, the lowest energy solutions are computed using the MPLP algorithm by Sontag et. al.
        First the energy .csv files are read and internally indexed. Then the required input files for the MPLP
        algorithm is written.
        The input files are then processed and the wanted number of solutions computed. The last step contains
        the creation of pdb files of the computed solutions, output .txt .html files and a final Pymol session
        containing all structures.

        Parameters
        ----------
        num_solutions: int
            Number of solutions to compute.
        ligand_scaling: int
            Scaling factor that will be multiplied with the ligand energies [default: 1]
        """
        from pocketoptimizer.preparation.aacodes import three2one
        from pocketoptimizer.solving import sontag_solver
        from pocketoptimizer.utility.energy_reader import EnergyReader
        from pocketoptimizer.utility.index_mapper import IndexMapper
        from pocketoptimizer.utility.sontag_writer import SontagWriter
        from pocketoptimizer.design.design_solution import DesignSolution
        from pocketoptimizer.design.html_reporter import HtmlReporter
        from pocketoptimizer.design.pymol_reporter import PymolReporter
        from pocketoptimizer.design.txt_reporter import TxtReporter

        if not self.mutations:
            logger.error('No mutations have been defined. Please define at least one mutation.')
            raise AttributeError('No mutations have been defined. Please define at least one mutation.')
        if not self.rotamer_path:
            logger.error('No rotamer library has been defined. Sample side chain rotamers first.')
            raise AttributeError('No rotamer library has been defined. Sample side chain rotamers first.')
        if not self.scorer:
            logger.error('No scoring method has been defined. Calculate energies first.')
            raise AttributeError('No scoring method has been defined. Calculate energies first.')

        mutations_copy = copy.deepcopy(self.mutations)

        sidechain_positions = {}
        for position in mutations_copy:
            pos = f'{position["chain"]}_{position["resid"]}'
            sidechain_positions[pos] = position["mutations"]

        design_mutations = ''
        for position, mutations in sidechain_positions.items():
            chain, resid = position.split('_')
            design_mutations += f'_{chain}{resid}'
            for mutation in mutations:
                design_mutations += f'{three2one[mutation]}'
        design_mutations = design_mutations[1:]
        if len(design_mutations) > 250:
            logger.warning('File name too long.')
            design_mutations = design_mutations[0:249]

        design_full_name = os.path.join(f'{self.forcefield}_{self.rotamer_path.split("/")[-1]}', design_mutations, f'{self.scorer}_scaling_{ligand_scaling}')
        solver_path = os.path.join(self.work_dir, 'solver', design_full_name)

        os.makedirs(solver_path, exist_ok=True)

        index_mapper = IndexMapper.from_conformer_files(rotamer_dir=self.rotamer_path, lig_poses=self.ligand_poses, mutations=self.mutations)
        index_mapper.set_scaling_factor('ligand', ligand_scaling)
        index_file_name = os.path.join(solver_path, 'index.dat')

        index_mapper.write_index_file(index_file_name)

        energy_reader = EnergyReader(work_dir=self.work_dir,
                                     mutations=mutations_copy,
                                     index_mapper=index_mapper,
                                     forcefield=self.forcefield,
                                     rotamer_path=self.rotamer_path,
                                     scorer=self.scorer)
        energy_reader.read_energies()

        sw = SontagWriter(index_mapper=index_mapper, energy_reader=energy_reader)
        sw.write_sontag(solver_path)

        os.makedirs(os.path.join(solver_path, 'solutions'), exist_ok=True)
        sontag_solver.calculate_design_solutions(
            solver_bin=self.settings.SOLVER_BIN,
            temp_dir=self.settings.TMP_DIR,
            out_path=solver_path,
            num_solutions=num_solutions,
            penalty_energy=1e10,
            exclude=None,
        )

        solution_file = os.path.join(solver_path, 'solutions', 'all_solutions.txt')
        index_file = os.path.join(solver_path, 'index.dat')
        energy_file = os.path.join(solver_path, 'lambdas.txt')

        logger.info('Parse calculated solutions.')
        solution = DesignSolution(work_dir=self.work_dir,
                                  solution_file=solution_file,
                                  index_file=index_file,
                                  energy_file=energy_file,
                                  mutations=mutations_copy,
                                  forcefield=self.forcefield,
                                  scaffold=self.built_scaffold,
                                  rotamer_path=self.rotamer_path,
                                  ligand=self.ligand_protonated,
                                  ligand_poses=self.ligand_poses,
                                  peptide=self.peptide,
                                  scorer=self.scorer)
        logger.info(f'Read {solution.get_solution_number()} solution(s) from solver output.')
        solution.read_detailed_self_energies()
        solution.read_detailed_pair_energies()

        design_outdir = os.path.join(self.work_dir, 'designs', design_full_name)
        os.makedirs(design_outdir, exist_ok=True)

        logger.info('Write txt report.')
        txt_reporter = TxtReporter(solution, sidechain_positions, design_outdir)
        txt_reporter.create_reports()
        logger.info('Wrote solution report text file(s).')
        txt_reporter.create_summary()
        logger.info('Wrote summary text file.')

        logger.info('Write html report.')
        html_reporter = HtmlReporter(design_solutions=solution, sidechain_positions=sidechain_positions,
                                     output_dir=design_outdir)
        html_reporter.create_reports()
        logger.info('Wrote solution report html file(s).')
        plot = True
        if solution.get_solution_number() < 1:
            plot = False
        if solution.get_solution_number() < 2:
            logo = False
        else:
            for position, residue in sidechain_positions.items():
                if solution.is_mutable(position):
                    logo = True
                    break
                else:
                    logo = False
        html_reporter.create_summary(plot, logo)
        logger.info('Wrote summary html file.')

        logger.info('Creating design structure files.')

        logger.info('Create Structures.')
        solution.create_design_structures(design_path=design_outdir)
        logger.info('Creating pymol scripts.')

        pymol_reporter = PymolReporter(solution, sidechain_positions, design_outdir)
        pymol_reporter.create_pymol_scripts(
            scaffold_pdb=self.built_scaffold,
            wt_ligand=self.ligand_protonated
        )

        logger.info('All solutions are processed.')

        logger.info(f'{num_solutions} best design solution(s) for design with forcefield: {self.forcefield}, scoring method: {self.scorer} and ligand scaling: {ligand_scaling} identified.')

    def design_multi(self, designs: List[Dict[str, Union[str, int]]]) -> NoReturn:
        """
        Wrapper function to run multiple designs.

        Parameters
        ----------
        designs: list
            Contains dictionaries containing the function arguments for the designs as keys, i.e.
            designs = [{'num_solutions':10, 'ligand_scaling':1]
        """
        for design in tqdm(designs, desc='Design'):
            self.design(num_solutions=design['num_solutions'],
                        ligand_scaling=design['ligand_scaling'])

    def clean(self, scaffold: bool = False, ligand: bool = False) -> None or NoReturn:
        """
        Remove all prepared scaffold and ligand files to start an entirely new design run.
        Delete settings.py, pocketoptimizer.log files and reset mutations, new DesignPipeline object needs to be initialized.

        Parameters
        ----------
        scaffold: bool
            Whether to remove prepared scaffold files [default: False]
        ligand: bool
            Whether to remove prepared ligand files [default: False]
        """
        import shutil

        def delete_prepared_scaffold_files(project_dir: str, ff: str) -> None or NoReturn:
            """
            Remove prepared scaffold files

            Parameters
            ----------
            project_dir: str
                Working directory
            ff: str
                forcefield used
            """
            dirs = ['scaffold', 'energies', 'solver', 'designs']
            paths = []
            if '..' in project_dir:
                logger.error('.. is not allowed in project_names.')
                return
            for dir in dirs:
                for file in glob.glob(os.path.join(project_dir, dir, '**', '*'), recursive=True):
                    if ff in file:
                        paths.append(file)
            if not len(paths):
                logger.info('No scaffold files found.')
            # Remove all files recursively in subdirectories
            for path in paths:
                if os.path.exists(path) and os.path.isfile(path):
                    os.remove(path)
            # Remove empty folders recursively
            for path in paths:
                if os.path.exists(path) and os.path.isdir(path):
                    shutil.rmtree(path)
            for dir in dirs:
                if os.path.exists(os.path.join(project_dir, dir)) and len(os.listdir(os.path.join(project_dir, dir))) == 0:
                    shutil.rmtree(os.path.join(project_dir, dir))

        def delete_prepared_ligand_files(project_dir: str, ff: str) -> None or NoReturn:
            """
            Remove prepared ligand files

            Parameters
            ----------
            project_dir: str
                Working directory
            ff: str
                forcefield used
            """
            paths = []
            if '..' in project_dir:
                logger.error('.. is not allowed in project_names.')
                return
            for file in glob.glob(os.path.join(project_dir, 'ligand', '**', '*'), recursive=True):
                if ff in file:
                    paths.append(file)
            if not len(paths):
                logger.info('No ligand files found.')
            # Remove all files recursively in subdirectories
            for path in paths:
                if os.path.exists(path) and os.path.isfile(path):
                    os.remove(path)
            # Remove empty folders recursively
            for path in paths:
                if os.path.exists(path) and os.path.isdir(path):
                    shutil.rmtree(path)
        if scaffold:
            delete_prepared_scaffold_files(project_dir=self.work_dir, ff=self.forcefield)
            logger.info('All scaffold files are deleted.')
        if ligand:
            delete_prepared_ligand_files(project_dir=self.work_dir, ff=self.forcefield)
            logger.info('All ligand files are deleted.')
        if os.path.isfile(os.path.join(self.work_dir, 'settings', 'settings.py')):
            os.remove(os.path.join(self.work_dir, 'settings', 'settings.py'))
            os.remove(os.path.join(self.work_dir, 'settings', '__init__.py'))
            rmtree(os.path.join(self.work_dir, 'settings'))
            logger.info('Deleted settings file.')
        if os.path.isfile(os.path.join(self.work_dir, os.environ.get('POCKETOPTIMIZER_LOGFILE'))):
            os.remove(os.path.join(self.work_dir, os.environ.get('POCKETOPTIMIZER_LOGFILE')))
            del os.environ['POCKETOPTIMIZER_LOGFILE']
            logger.info('Deleted log file.')
        # Reset self.mutations, since new design run will be initialized
        self.mutations = []
        logger.info('All files were deleted.')
        logger.info('Initialize a new DesignPipeline.')


# CLI
if __name__ == '__main__':
    import sys
    import argparse
    sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

    import pathlib
    # Get path of working directory
    working_dir = str(pathlib.Path().resolve())

    import pocketoptimizer
    parser = argparse.ArgumentParser(description='PocketOptimizer computational protein design pipeline CLI, for more options use API.')
    parser.add_argument('-ff', '--forcefield', type=str, nargs=1, help='Force field to be used either: amber_ff14SB or charmm36', default=['amber_ff14SB'], required=False)
    parser.add_argument('-r', '--receptor', type=str, nargs=1, help='Protein input structure file in pdb format', required=True)
    parser.add_argument('-l', '--ligand', type=str, nargs=1, help='Ligand input structure file', required=True)
    parser.add_argument('--ph', type=float, nargs=1, default=[7.0], help='ph value for side chain and ligand protonation', required=False)
    parser.add_argument('--keep_chains', type=str, nargs='*', help='Chains to keep by their chain identifiers', required=False)
    parser.add_argument('--discard_mols', type=str, nargs='*', help='Special molecules to exclude by their chain and residue identifier (A:1), '
                                                                           'per default everything, but peptides have to be defined manually', required=False)
    parser.add_argument('--mutations', type=str, nargs='+', help='Mutations (A:1:ALA)', required=True)
    parser.add_argument('--peptide_mutations', type=str, nargs='*', default=None, help='Peptide mutations (1:ALA)', required=False)
    parser.add_argument('--flex_peptide_res', type=str, nargs='*', default=None, help='Peptide residues to sample flexiblity', required=False)
    parser.add_argument('--vdw_thresh', type=float, nargs=1, default=[100.0], help='VdW-energy threshold for rotamer and ligand pose sampling (kcal/mol)', required=False)
    parser.add_argument('--rot_lib', type=str, nargs=1, default=['dunbrack'], help='Rotamer library, options are: dunbrack or cmlib', required=False)
    parser.add_argument('--nconfs', type=int, nargs=1, default=[50], help='Number of ligand conformers to sample', required=False)
    parser.add_argument('--lig_rot', '--lig_rot', type=float, nargs=1, default=[20], help='Maximum ligand rotation', required=False)
    parser.add_argument('--lig_rot_steps', '--lig_rot_steps', type=float, nargs=1, default=[20], help='Ligand rotation steps', required=False)
    parser.add_argument('--lig_trans', '--lig_trans', type=float, nargs=1, default=[1], help='Maximum ligand translation', required=False)
    parser.add_argument('--lig_trans_steps', '--lig_trans_steps', type=float, nargs=1, default=[0.5], help='Ligand translation steps', required=False)
    parser.add_argument('--max_poses', '--max_poses', type=int, nargs=1, default=[10000], help='Maximum number of ligand poses to sample', required=False)
    parser.add_argument('--sampling_pocket', type=str, nargs=1, default=['GLY'], help='Sampling pocket for rotamer and ligand pose sampling', required=False)
    parser.add_argument('--scoring', type=str, nargs=1, default=['vina'], help='Scoring function, options are: vina, vinardo, dkoes_scoring, ad4_scoring, force_field', required=False)
    parser.add_argument('--lig_scaling', type=int, nargs=1, default=[1], help='Ligand scaling factor', required=False)
    parser.add_argument('--num_solutions', type=int, nargs=1, default=[10], help='Number of design solutions to calculate', required=False)
    parser.add_argument('--ncpus', type=int, nargs=1, default=[1], help='Number of CPUs for multiprocesing', required=False)
    parser.add_argument('--cuda', type=int, nargs=1, default=[0], help='Enabling cuda for GPU based minimization', required=False)

    # cuda, conf_method, vdw

    # Custom class for cleaning the working directory
    class CleanWorkDir(argparse.Action):
        def __call__(self, parser, namespace, values, option_string):
            design = pocketoptimizer.DesignPipeline(work_dir=working_dir, forcefield=values[0])
            design.clean(scaffold=True, ligand=True)
            parser.exit() # exits the program with no more arg parsing and checking

    parser.add_argument('--clean', nargs=1, action=CleanWorkDir, help='Clean the working directory', required=False)
    args = parser.parse_args()

    discard_mols = []
    if args.discard_mols is not None:
        for mol in args.discard_mols:
            try:
                chain, resid = mol.split(':')
                discard_mols.append({'chain': chain, 'resid': resid})
            except ValueError:
                logger.error('Define molecules in the following format CHAIN:RESID')
                raise argparse.ArgumentTypeError('Define molecules in the following format CHAIN:RESID')

    mutations = []
    for mutation in args.mutations:
        try:
            chain, resid, resname = mutation.split(':')
            mutations.append({'chain': chain, 'resid': resid, 'mutations': [resname]})
        except ValueError:
            logger.error('Define mutations in the following format CHAIN:RESID:RESNAME')
            raise argparse.ArgumentTypeError('Define mutations in the following format CHAIN:RESID:RESNAME')

    if args.peptide_mutations:
        peptide_mutations = []
        for mutation in args.peptide_mutations:
            try:
                resid, resname = mutation.split(':')
                peptide_mutations.append({'resid': resid, 'mutation': resname})
            except ValueError:
                logger.error('Define mutations in the following format RESID:RESNAME')
                raise argparse.ArgumentTypeError('Define mutations in the following format RESID:RESNAME')

    # Initialize new DesignPipeline in current working directory.
    if not args.flex_peptide_res:
        design = pocketoptimizer.DesignPipeline(work_dir=working_dir, forcefield=args.forcefield[0], ph=args.ph[0], ncpus=args.ncpus[0], peptide=False)
        design.parameterize_ligand(input_ligand=args.ligand[0])
        design.prepare_lig_conformers(nconfs=args.nconfs[0])
        design.prepare_protein(protein_structure=args.receptor[0], keep_chains=args.keep_chains, backbone_restraint=True,
                               cuda=bool(args.cuda[0]), discard_mols=discard_mols)
    else:
        design = pocketoptimizer.DesignPipeline(work_dir=working_dir, forcefield=args.forcefield[0], ph=args.ph[0], ncpus=args.ncpus[0], peptide=True)
        design.prepare_protein(protein_structure=args.receptor[0], keep_chains=args.keep_chains, backbone_restraint=True,
                               discard_mols=discard_mols, peptide_structure=args.ligand[0], peptide_mutations=peptide_mutations, cuda=bool(args.cuda[0]))
        design.prepare_peptide_conformers(positions=[resid for resid in args.flex_peptide_res], vdw_filter_thresh=args.vdw_thresh[0])

    design.set_mutations(mutations)
    design.prepare_mutants(sampling_pocket=args.sampling_pocket[0])
    design.sample_sidechain_rotamers(library=args.rot_lib[0], vdw_filter_thresh=args.vdw_thresh[0])
    design.sample_lig_poses(method='grid', grid={'trans': [args.lig_trans[0], args.lig_trans_steps[0]], 'rot': [args.lig_rot[0], args.lig_rot_steps[0]]} ,
                            vdw_filter_thresh=args.vdw_thresh[0], max_poses=args.max_poses[0])
    design.calculate_energies(scoring=args.scoring[0])
    design.design(num_solutions=args.num_solutions[0], ligand_scaling=args.lig_scaling[0])
