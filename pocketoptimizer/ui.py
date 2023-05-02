import copy
import shutil

from tqdm.auto import tqdm
import os
import glob
from shutil import which
from typing import List, Dict, Union, NoReturn
import logging.config

logger = logging.getLogger(__name__)

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

    def __init__(self, work_dir: str, forcefield: str, ph: float = 7.0, ncpus: int = 1):
        """
        Constructor method

        Parameters
        ----------
        work_dir: str
            Working directory, all files will be created in this directory
        forcefield: str
            Force field used for sampling procedures and energy calculation steps [options: 'amber_ff14SB', 'charmm36']
        elec: float
            Scaling factor for electrostatic components [default: 0.01]
        pH: float
            pH value that will be used to adjust the protonation states [default: 7.0]
        ncpus:
            Number of CPUs to use [default: 1]
        """
        import tempfile
        import multiprocessing as mp
        from pocketoptimizer.path import path

        self.work_dir = work_dir
        os.chdir(self.work_dir)

        self.logfile = os.path.join(self.work_dir, 'pocketoptimizer.log')
        # Set environment variable to standard path
        os.environ["POCKETOPTIMIZER_LOGFILE"] = self.logfile
        # Set environment variables for MATCH
        os.environ["MATCH"] = os.path.abspath(os.path.join(path(), '..', 'MATCH_RELEASE', 'MATCH'))
        os.environ["PerlChemistry"] = os.path.abspath(os.path.join(path(), '..', 'MATCH_RELEASE', 'PerlChemistry'))

        logging.config.fileConfig(
            os.path.join(path(), "logging.ini"), disable_existing_loggers=False
        )
        logger.info(f'Logging to: {self.logfile}')

        if not os.path.isdir(self.work_dir):
            logger.error(f'{self.work_dir} does not exist.')
            raise FileNotFoundError(f'{self.work_dir} does not exist.')

        self.tmp_dir = tempfile.gettempdir()

        self.match = os.path.abspath(os.path.join(path(), '..', 'MATCH_RELEASE', 'MATCH', 'scripts', 'MATCH.pl'))
        self.solver = path(bin_file='sontag_solver')
        self.psfgen = path(bin_file='psfgen')

        obabel_path = which('obabel')
        if not obabel_path:
            logger.warning("Obabel not found make sure it is installed and in your path.")
        self.obabel = obabel_path

        antechamber_path = which('antechamber')
        if not antechamber_path:
            antechamber_path = '#'
            logger.warning("Antechamber not found make sure ambertools is installed and in your path.")
        self.antechamber = antechamber_path
        tleap_path = which('tleap')
        if not tleap_path:
            tleap_path = '#'
            logger.warning("Tleap not found make sure ambertools is installed and in your path.")
        self.tleap = tleap_path
        parmchk2_path = which('parmchk2')
        if not parmchk2_path:
            parmchk2_path = '#'
            logger.warning("Parmchk2 not found make sure it is installed and in your path.")
        self.parmchk2 = parmchk2_path
        smina_path = which('smina')
        if not smina_path:
            smina_path = '#'
            logger.warning("Smina not found make sure it is installed and in your path.")
        self.smina = smina_path

        self.built_scaffold = ''
        self.prepared_protein = ''
        self.ligand_protonated = ''
        self.built_ligand_params = {}
        self.ligand_conformers = ''
        self.ligand_poses_pdb = ''
        self.ligand_poses_xtc = ''
        self.mutations = []
        self.scoring = {
            'smina': ['vina', 'vinardo', 'ad4_scoring'],
            'ff': ['amber_ff14SB', 'charmm36'],
        }
        self.scorer = ''
        self.forcefield = ''
        self._set_ff(forcefield)
        self.library = ''
        self.rotamer_path = ''
        self.ph = ph
        self.temperature = 300.0
        self.ncpus = ncpus
        if self.ncpus > mp.cpu_count():
            logger.warning(f'More CPUs defined than available, setting maximum of {mp.cpu_count()} CPUs.')
            self.ncpus = mp.cpu_count()

    def _set_ff(self, forcefield: str) -> NoReturn:
        """
        Sets protein and ligand force fields used for all sampling procedures and energy calculations.

        Parameters
        ----------
        forcefield: str
            Force field. Options are 'amber_ff14SB' and 'charmm36' (Default: 'amber_ff14SB').
         """
        if forcefield in self.scoring['ff']:
            self.forcefield = forcefield
            self.built_scaffold = os.path.join(self.work_dir, 'scaffold', self.forcefield, 'scaffold.pdb')
            self.prepared_protein = os.path.join(self.work_dir, 'scaffold', self.forcefield, 'protein_preparation', 'prepared_scaffold.pdb')
            self.ligand_protonated = os.path.join(self.work_dir, 'ligand', self.forcefield, 'ligand.mol2')
            self.built_ligand_params = {'params_folder': os.path.join(self.work_dir, 'ligand', self.forcefield, 'params'),
                                        'mol2': [os.path.join(self.work_dir, 'ligand', self.forcefield, 'params','ligand.mol2')],
                                        'frcmod': [os.path.join(self.work_dir, 'ligand', self.forcefield, 'params','ligand.frcmod')],
                                        'prm': [os.path.join(self.work_dir, 'ligand', self.forcefield, 'params','ligand.prm')],
                                        'rtf': [os.path.join(self.work_dir, 'ligand', self.forcefield, 'params','ligand.rtf')]}
            self.ligand_conformers = os.path.join(self.work_dir, 'ligand', self.forcefield, 'conformers', 'ligand_confs.pdb')
            self.ligand_poses_pdb = os.path.join(self.work_dir, 'ligand', self.forcefield, 'poses', 'ligand_poses.pdb')
            self.ligand_poses_xtc = os.path.join(self.work_dir, 'ligand', self.forcefield, 'poses', 'ligand_poses.xtc')
        else:
            logger.error(f'{forcefield} not supported. Use {self.scoring["ff"][0]} or {self.scoring["ff"][1]}.')
            raise ValueError(f'{forcefield} not supported. Use {self.scoring["ff"][0]} or {self.scoring["ff"][1]}.')

    def _set_rot_lib(self, library: str) -> NoReturn:
        """
        Sets the rotamer path, based on the rotamer libary

        Parameters
        ----------
        library: str
            Rotamer library used for rotamer sampling.
            Options are either: 'dunbrack' or 'cmlib'
        """
        if library not in ['dunbrack', 'cmlib']:
            logger.error('Rotamer library not implemented.')
            raise NotImplementedError('Rotamer library not implemented.')
        self.library = library
        self.rotamer_path = os.path.join(self.work_dir, 'scaffold', self.forcefield, 'rotamers', self.library)

    def _set_temperature(self, temperature: float) -> NoReturn:
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
                                               mutations=mutations,
                                               forcefield=self.forcefield)
        self.mutations = mutation_processor.process_mutations()

    def parameterize_ligand(self, input_ligand: str, addHs: bool = True) -> NoReturn:
        """
        Interface for Ligand parameterization.
        Adds hydrogen atoms to the input molecule, produces GAFF2 or CGenFF_2b6 parameters.

        Parameters
        ----------
        input_ligand: str
            Path to a ligand file
        addHs: bool
            Whether to add hydrogen atoms [default: True]
        """
        from pocketoptimizer.preparation.structure_building import SystemBuilder

        if not os.path.isfile(self.ligand_protonated):

            system = SystemBuilder(design_pipeline=self,
                                   structure=os.path.join(self.work_dir, input_ligand))

            system.parameterize_ligand(addHs=addHs)

        else:
            logger.info('Ligand is already parametrized.')

    def prepare_protein(self, protein_structure: str, keep_chains: List = [], discard_mols: List[Dict[str, str]] = [], backbone_restraint: bool = True, cuda: bool = False) -> NoReturn:
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
            discard_mols: list
                List of dictionaries containing chain and resid of the molecules to discard. In the case of amino acid ligands,
                these have to be manually discarded with this option
            backbone_restraint: bool
                Restraints the backbone during minimization. [default: True]
            cuda: bool
                If the minimization should be performed on a GPU. [default: False]
            """
            from pocketoptimizer.preparation.structure_building import SystemBuilder

            logger.info('Start Protein Preparation.')

            if not type(keep_chains) == list:
                logger.error('Define chains in a list.')
                raise TypeError('Define chains in a list.')

            if not type(discard_mols) == list:
                logger.error('Define molecules in a list of dictionaries.')
                raise TypeError('Define molecules in a list of dictionaries.')

            system = SystemBuilder(design_pipeline=self,
                                   structure=os.path.join(self.work_dir, protein_structure))

            os.makedirs(os.path.dirname(self.prepared_protein), exist_ok=True)
            system.prepare_protein(keep_chains=keep_chains,
                                   discard_mols=discard_mols)

            # Check if there is a built scaffold
            if not os.path.isfile(self.built_scaffold):
                # Minimize the structure with ligand in pocket
                logger.info('Building complex.')
                # Build the complex before minimization
                system.build_complex()
                input_path = os.path.join(self.work_dir, 'scaffold', self.forcefield, 'protein_params', 'native_complex')

                system.minimize_structure(
                    structure_path=input_path,
                    cuda=cuda,
                    restraint_bb=backbone_restraint)

                logger.info('Your protein was successfully minimized and can be used for design now.')

            else:
                logger.info('Your scaffold is already prepared.')

            logger.info('Protein preparation finished.')

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
            RMSD cutoff, only used together with the 'confab' method. [default: 0.5 Angstrom]
        ecutoff: float
            Energy cutoff, only used together with the 'confab' method. [default: 50 kcal/mol]
        """
        from pocketoptimizer.sampling import conformer_generator_obabel

        os.makedirs(os.path.join(self.work_dir, 'ligand', self.forcefield, 'conformers'), exist_ok=True)

        if method == 'genetic' or method == 'confab':
            conformer_generator_obabel.conformer_generator(
                obabel_path=self.obabel, infile=self.ligand_protonated, conf_file_name=self.ligand_conformers,
                method=method, nconfs=nconfs, score=score, rcutoff=rcutoff, ecutoff=ecutoff)
        else:
            logger.error('Conformer generation method not defined! Try genetic or confab.')
            raise ValueError('Conformer generation method not defined! Try genetic or confab.')

    def prepare_mutants(self, sampling_pocket: str = 'ALA') -> NoReturn:
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

        if not self.mutations:
            logger.error('No mutations have been defined. Please define at least one mutation.')
            raise ValueError('No mutations have been defined. Please define at least one mutation.')

        logger.info('Start building mutated protein scaffold variants.')

        if sampling_pocket not in aa:
            logger.error('Amino acid unknown.')
            raise ValueError('Amino acid unknown.')
        logger.info(f'Build {sampling_pocket} sampling pockets.')

        system = SystemBuilder(design_pipeline=self, structure=self.built_scaffold)
        # Build all single and pairwise mutated scaffolds and all sampling pockets for sampling rotamers and ligand poses
        system.build_complex(sampling_pocket=sampling_pocket)
        logger.info('Scaffold building done.')

    def sample_sidechain_rotamers(self, library: str = 'dunbrack', vdw_filter_thresh: float = 100.0,
                                  dunbrack_filter_thresh: float = 0.01, expand: List[str] = ['chi1', 'chi2'],
                                  accurate: bool = False, include_native: bool = True) -> NoReturn:
        """Uses the dunbrack or cm_lib rotamer library to compute possible rotamers of defined residues.
        Energies for rotamer pruning are calculated using ffevaluate with the amber_ff14SB or charmm36 force field


        Parameters
        ----------
        library: str
            Rotamer library, options are: 'cm_lib' or 'dunbrack', 'dunbrack' is recommended [default: 'dunbrack']
        vdw_filter_thresh: float
            Threshold value used to filter rotamers [default: 100.0 kcal/mol]
        dunbrack_filter_thresh: float
            Filter threshold, rotamers having probability of occurence lower than filter threshold will
            be pruned if their rotameric mod does occur more than once, value should be between 0 and 1
            (-1: no pruning, 1: pruning of all rotamers with duplicate rotamer modes) [default: 1 %]
        expand: list
            List of which chi angles to expand, [default: ['chi1', 'chi2']]
        accurate: bool
            Whether to expand chi-angles by +/-1 std or also +/-0.5 std [default: False]
        include_native: bool
            Whether to include the native rotamers from the structure [default: True]
        """

        from pocketoptimizer.sampling.sidechain_rotamers_ffev import FFRotamerSampler

        if not self.mutations:
            logger.error('No mutations have been defined. Please define at least one mutation.')
            raise ValueError('No mutations have been defined. Please define at least one mutation.')
        # set rotamer path for energy calculations
        self._set_rot_lib(library)

        if vdw_filter_thresh < 0:
            vdw_filter_thresh = 0

        rotamer_sampler = FFRotamerSampler(design_pipeline=self)

        rotamer_sampler.rotamer_sampling(
            vdw_filter_thresh=vdw_filter_thresh,
            dunbrack_prob=dunbrack_filter_thresh,
            expand=expand,
            accurate=accurate,
            include_native=include_native)
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
        vdw_filter_thresh: float
            Energy threshold for filtering ligand poses [default: 100 kcal/mol]
        max_poses: int
            Maximum amount of poses to produce [default: 10000].
        """
        from pocketoptimizer.sampling import ligand_poses

        if method == 'grid' and not grid:
            logger.error('Define grid, if grid method should be used.')
            raise ValueError('Define grid, if grid method should be used.')

        if vdw_filter_thresh < 0:
            vdw_filter_thresh = 0

        sampling_scaffold_ligand = os.path.join(self.work_dir, 'scaffold', self.forcefield, 'protein_params',
                                                'ligand_sampling_pocket', 'structure.pdb')

        if os.path.isfile(sampling_scaffold_ligand):
            if not os.path.isfile(self.ligand_poses_pdb):
                logger.info('Sample ligand poses.')
                sampler = ligand_poses.PoseSampler(design_pipeline=self,
                                                   max_poses=max_poses,
                                                   sampling_method=method)

                sampler.create_ligand_ensemble(grid=grid,
                                               vdw_filter_thresh=vdw_filter_thresh)
            else:
                logger.info('Ligand poses are already sampled.')
        else:
            logger.error('Ligand pose sampling scaffold has not been prepared yet.')

    def calculate_energies(self, scoring: str = 'vina') -> NoReturn:
        """Computes self and pairwise energies of all ligand and scaffold components.

        The energy calculation procedure consists of three parts.
        1. Sidechain vs. Scaffold energies
        2. Sidechain vs.Sidechain
        3. Ligand vs. Scaffold
        4. Ligand vs. Sidechain

        The former two are related to force field computations while the latter two can be carried out using a
        specific scoring function.
        Output energy scores will be written into .csv files under /work_dir/energies/ directory.

        Parameters
        ----------
        scoring: str
            Method that will be used for scoring [default: 'vina']
            The following methods are available: vina, vinardo, ad4_scoring, amber_ff14SB, charmm36
        """
        from pocketoptimizer.scoring.sidechain_scaffold_energies import SidechainSelfScorer
        from pocketoptimizer.scoring.sidechain_pair_energies import SidechainPairScorer

        os.chdir(self.work_dir)

        if not self.mutations:
            logger.error('No mutations have been defined. Please define at least one mutation.')
            raise ValueError('No mutations have been defined. Please define at least one mutation.')
        if not self.library:
            logger.error('No rotamer library has been defined. Sample side chain rotamers first.')
            raise ValueError('No rotamer library has been defined. Sample side chain rotamers first.')

        if scoring in self.scoring['ff']:
            if scoring != self.forcefield:
                logger.warning(f'Can not score with {scoring} if force field is {self.forcefield}.')
                raise RuntimeWarning(f'Can not score with {scoring} if force field is {self.forcefield}.')
        if scoring not in self.scoring['smina'] and scoring not in self.scoring['ff']:
            logger.error(f'Scorer: {scoring} is not implemented.')
            raise NotImplementedError(f'Scorer: {scoring} is not implemented.')
        self.scorer = scoring

        logger.info('Start energy calculations.')
        logger.info(f'Using {self.ncpus} CPUs for multiprocessing.')

        logger.info('Calculate Sidechain-Scaffold Energies.')
        self_scorer = SidechainSelfScorer(design_pipeline=self)
        self_scorer.calculate_scaffold()

        logger.info('Calculate Sidechain-Pair Energies.')
        pair_scorer = SidechainPairScorer(design_pipeline=self)
        pair_scorer.calculate_pairs()

        logger.info('Calculate Ligand-Scaffold/Sidechain-Interaction-Energies.')
        if self.scorer in self.scoring['smina']:
            from pocketoptimizer.scoring.smina_scorer import SminaScorer
            smina_scorer = SminaScorer(design_pipeline=self)

            smina_scorer.run_smina_scorer()
        elif self.scorer in self.scoring['ff']:
            from pocketoptimizer.scoring.ff_scorer import LigandScorer

            ff_scorer = LigandScorer(design_pipeline=self)
            ff_scorer.run_ff_scorer()

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
        algorithm are written.
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
            raise ValueError('No mutations have been defined. Please define at least one mutation.')
        if not self.library:
            logger.error('No rotamer library has been defined. Sample side chain rotamers first.')
            raise ValueError('No rotamer library has been defined. Sample side chain rotamers first.')
        if not self.scorer:
            logger.error('No scoring method has been defined. Calculate energies first.')
            raise ValueError('No scoring method has been defined. Calculate energies first.')
        if num_solutions < 1:
            logger.error('Need at least one solution.')
            raise ValueError('Need at least one solution.')

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

        design_full_name = os.path.join(f'{self.forcefield}_{self.library}', design_mutations, f'{self.scorer}_scaling_{ligand_scaling}')
        solver_path = os.path.join(self.work_dir, 'solver', design_full_name)

        os.makedirs(solver_path, exist_ok=True)

        index_mapper = IndexMapper.from_conformer_files(design_pipeline=self)
        index_mapper.set_scaling_factor('ligand', ligand_scaling)
        index_file_name = os.path.join(solver_path, 'index.dat')

        index_mapper.write_index_file(index_file_name)

        energy_reader = EnergyReader(design_pipeline=self,
                                     index_mapper=index_mapper)
        energy_reader.read_energies()

        sw = SontagWriter(index_mapper=index_mapper, energy_reader=energy_reader)
        sw.write_sontag(solver_path)

        os.makedirs(os.path.join(solver_path, 'solutions'), exist_ok=True)
        sontag_solver.calculate_design_solutions(solver_bin=self.solver,
                                                 design_pipeline=self,
                                                 temp_dir=self.tmp_dir,
                                                 out_path=solver_path,
                                                 num_solutions=num_solutions,
                                                 penalty_energy=1e10,
                                                 exclude=None)

        solution_file = os.path.join(solver_path, 'solutions', 'all_solutions.txt')
        index_file = os.path.join(solver_path, 'index.dat')
        energy_file = os.path.join(solver_path, 'lambdas.txt')

        logger.info('Parse calculated solutions.')
        solution = DesignSolution(design_pipeline=self,
                                  solution_file=solution_file,
                                  index_file=index_file,
                                  energy_file=energy_file)
        logger.info(f'Read {solution.get_solution_number()} solution(s) from solver output.')
        solution.read_detailed_self_energies()
        solution.read_detailed_pair_energies()
        # Remove solver directory
        shutil.rmtree(os.path.join(self.work_dir, 'solver'))

        design_outdir = os.path.join(self.work_dir, 'designs', design_full_name)
        os.makedirs(design_outdir, exist_ok=True)

        logger.info('Write text report.')
        txt_reporter = TxtReporter(solution, sidechain_positions, design_outdir)
        txt_reporter.create_reports()
        logger.info('Wrote solution report text file(s).')
        txt_reporter.create_summary()
        logger.info('Wrote summary text file.')

        logger.info('Write html report.')
        html_reporter = HtmlReporter(design_solutions=solution,
                                     sidechain_positions=sidechain_positions,
                                     output_dir=design_outdir)
        html_reporter.create_reports()
        logger.info('Wrote solution report html file(s).')
        html_reporter.create_summary()
        logger.info('Wrote summary html file.')

        logger.info('Creating design structure files.')

        logger.info('Create Structures.')
        solution.create_design_structures(design_path=design_outdir)
        logger.info('Creating PyMol scripts.')

        pymol_reporter = PymolReporter(solution, sidechain_positions, design_outdir)
        pymol_reporter.create_pymol_scripts(
            scaffold_pdb=self.built_scaffold,
            wt_ligand=self.ligand_protonated
        )

        logger.info(f'{solution.get_solution_number()} best design solution(s) for design with forcefield: {self.forcefield}, scoring method: {self.scorer} and ligand scaling: {ligand_scaling} identified.')

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
        Delete logfile and reset mutations, new DesignPipeline object needs to be initialized.

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
            logger.info('All scaffold files were deleted.')
        if ligand:
            delete_prepared_ligand_files(project_dir=self.work_dir, ff=self.forcefield)
            logger.info('All ligand files were deleted.')
        if os.path.isfile(self.logfile):
            os.remove(self.logfile)
            del os.environ['POCKETOPTIMIZER_LOGFILE']
            logger.info('Logfile was deleted.')


def main():
    import sys
    import argparse
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

    import pathlib
    # Get path of working directory
    working_dir = str(pathlib.Path().resolve())

    import pocketoptimizer
    parser = argparse.ArgumentParser(description='PocketOptimizer computational protein design pipeline CLI, for more options use API.')
    parser.add_argument('-ff', '--forcefield', type=str, help='Force field to be used either: amber_ff14SB or charmm36', default='amber_ff14SB', required=False)
    parser.add_argument('-r', '--receptor', type=str, help='Protein input structure file in pdb format', required=True)
    parser.add_argument('-l', '--ligand', type=str, help='Ligand input structure file', required=True)
    parser.add_argument('--ph', type=float, default=7.0, help='ph value for side chain and ligand protonation', required=False)
    parser.add_argument('--keep_chains', type=str, nargs='*', help='Chains to keep by their chain identifiers', required=False)
    parser.add_argument('--min_bb', action='store_true', help='Whether to minimize the proteins backbone')
    parser.add_argument('--discard_mols', type=str, nargs='*', help='Special molecules to exclude by their chain and residue identifier (A:1), '
                                                                           'per default everything, but peptides have to be defined manually', required=False)
    parser.add_argument('--mutations', type=str, nargs='+', help='Mutations (A:1:ALA)', required=True)
    parser.add_argument('--vdw_thresh', type=float, default=100.0, help='Energy threshold for rotamer and ligand pose sampling (kcal/mol)', required=False)
    parser.add_argument('--library', type=str, default='dunbrack', help='Rotamer library, options are: dunbrack or cmlib', required=False)
    parser.add_argument('--dunbrack_filter_thresh', type=float, default=0.01, help='Filter threshold for dunbrack rotamer library (between 0 and 1), default: 0.01', required=False)
    parser.add_argument('--accurate', action='store_true', help='Sample additional rotamers')
    parser.add_argument('--include_native', action='store_true', help='Include native rotamer')
    parser.add_argument('--nconfs', type=int, default=50, help='Number of ligand conformers to sample, default: 50', required=False)
    parser.add_argument('--rot', '--rot', type=float, default=20, help='Maximum ligand rotation, default: 20°', required=False)
    parser.add_argument('--rot_steps', '--rot_steps', type=float, default=20, help='Ligand rotation steps, default: 20°', required=False)
    parser.add_argument('--trans', '--trans', type=float, default=1, help='Maximum ligand translation, default: 1 Å', required=False)
    parser.add_argument('--trans_steps', '--trans_steps', type=float, default=0.5, help='Ligand translation steps, default 0.5 Å', required=False)
    parser.add_argument('--max_poses', '--max_poses', type=int, default=10000, help='Maximum number of ligand poses to sample, default: 10000', required=False)
    parser.add_argument('--sampling_pocket', type=str, default='ALA', help='Sampling pocket for rotamer and ligand pose sampling, default: ALA', required=False)
    parser.add_argument('--scoring', type=str, default='vina', help='Scoring function, options are: vina, vinardo, ad4_scoring, amber_ff14SB or charmm36', required=False)
    parser.add_argument('--scaling', type=int, default=1, help='Ligand scaling factor, default: 1', required=False)
    parser.add_argument('--num_solutions', type=int, default=10, help='Number of design solutions to calculate, default 10', required=False)
    parser.add_argument('--ncpus', type=int, default=1, help='Number of CPUs for multiprocessing', required=False)
    parser.add_argument('--cuda', action='store_true', help='Enabling cuda for GPU-based minimization')

    # Custom class for cleaning the working directory
    class CleanWorkDir(argparse.Action):
        def __call__(self, parser, namespace, values, option_string):
            design = pocketoptimizer.DesignPipeline(work_dir=working_dir, forcefield=values[0])
            design.clean(scaffold=True, ligand=True)
            parser.exit() # exits the program with no more arg parsing and checking

    parser.add_argument('--clean', nargs=1, action=CleanWorkDir, help='Clean the working directory', required=False)
    args = parser.parse_args()

    keep_chains = []
    if args.keep_chains:
        for chain in args.keep_chains:
            keep_chains.append(chain)

    discard_mols = []
    if args.discard_mols:
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


    # Initialize new DesignPipeline in current working directory
    design = pocketoptimizer.DesignPipeline(work_dir=working_dir,
                                            forcefield=args.forcefield,
                                            ph=args.ph,
                                            ncpus=args.ncpus)
    design.parameterize_ligand(input_ligand=args.ligand)
    design.prepare_lig_conformers(nconfs=args.nconfs)
    design.prepare_protein(protein_structure=args.receptor,
                           keep_chains=keep_chains,
                           backbone_restraint=not args.min_bb,
                           discard_mols=discard_mols,
                           cuda=args.cuda)

    design.set_mutations(mutations)
    design.prepare_mutants(sampling_pocket=args.sampling_pocket)
    design.sample_sidechain_rotamers(library=args.library,
                                     vdw_filter_thresh=args.vdw_thresh,
                                     dunbrack_filter_thresh=args.dunbrack_filter_thresh,
                                     accurate=args.accurate,
                                     include_native=args.include_native)
    design.sample_lig_poses(method='grid',
                            grid={'trans': [args.trans, args.trans_steps], 'rot': [args.rot, args.rot_steps]},
                            vdw_filter_thresh=args.vdw_thresh,
                            max_poses=args.max_poses)
    design.calculate_energies(scoring=args.scoring)
    design.design(num_solutions=args.num_solutions,
                  ligand_scaling=args.scaling)


if __name__ == '__main__':
    main()
