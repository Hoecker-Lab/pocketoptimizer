import itertools
import os
import shutil
import subprocess
from typing import List, Dict, Union, NoReturn
import logging
import htmd.builder.amber as amber
import htmd.builder.charmm as charmm
import numpy as np
from moleculekit.molecule import Molecule
from moleculekit.projections.metricrmsd import MetricRmsd
from moleculekit.tools.autosegment import autoSegment2
from moleculekit.tools.preparation import systemPrepare
from openmm import app
import openmm as mm
from simtk import unit
import parmed

from pocketoptimizer.utility.utils import fix_parameters, load_ff_parameters, create_pairs, Storer
from pocketoptimizer.preparation.aacodes import special2standard
from pocketoptimizer.path import path

logger = logging.getLogger(__name__)


class SystemBuilder(Storer):
    """
    Class providing system building functionalities to create systems for force field energy computations and minimization.
    """
    def __init__(self, structure: str, **kwargs):
        """
        Constructor method.

        Parameters
        ----------
        structure: str
            Path to structure
        """
        super().__init__(**kwargs)
        self.structure = structure

    def prepare_protein(self, keep_chains: List[str], discard_mols: List[Dict[str, str]] = []) -> NoReturn:
        """
        Protein cleaning and protonation procedure using Moleculekit.

        The wanted protein chain gets extracted and cleaned from ions, molecules etc.
        Gaps in the protein are closed and the protein is protonated

        Parameters
        ----------
        keep_chains: str
            Protein chain that will be extracted.
        discard_mols: list
            List of dictionaries containing chain and resid of the molecules to discard.
            In the case of amino acid ligands, these have to be manually discarded with this option. [default: []]
        """
        # check if protein is already prepared
        if not os.path.isfile(self.prepared_protein):
            renumber = False
            with open(self.structure, 'r') as origin_scaffold:
                for line in origin_scaffold:
                    if line.startswith('ATOM'):
                        try:
                            int(line[22:26])
                        except ValueError:
                            renumber = True
                            break

            prot = Molecule(self.structure)

            if renumber:
                self.renumber_residues(prot=prot)

            # Filter ligands, ions, waters out
            filter_str = f'protein'
            if discard_mols:
                for mol in discard_mols:
                        filter_str += f' and not (chain {mol["chain"]} and resid {mol["resid"]})'

            if keep_chains:
                selection_str = '('
                for i, chain in enumerate(keep_chains):
                    if i == len(keep_chains) - 1:
                        selection_str += f'chain {chain})'
                    else:
                        selection_str += f'chain {chain} or '
                filter_str = selection_str + ' and ' + filter_str

            prot.filter(filter_str, _logger=False)

            # Adding hydrogens and protonating
            logger.info(f'Protonate protein according to pH: {self.ph}.')
            prot_prep, prep_data = systemPrepare(prot, pH=self.ph, return_details=True, _molkit_ff=False)

            os.makedirs(os.path.split(self.prepared_protein)[0], exist_ok=True)
            # Finding gaps so gaps can get caps
            prot_prep = autoSegment2(prot_prep, _logger=False)

            prot_prep.write(self.prepared_protein, type='pdb')
            prep_data.to_excel(os.path.join(os.path.split(self.prepared_protein)[0], 'report.xlsx'))

            if not os.path.isfile(self.prepared_protein):
                logger.error('Protein preparation failed.')
                raise FileNotFoundError('Protein preparation failed.')

            # change scaffold name to prepared protein since new scaffold was prepared for building
            self.structure = self.prepared_protein

            logger.info('Successfully prepared protein structure.')
        else:
            logger.info('Protein is already prepared.')

    def renumber_residues(self, prot: Molecule):
        """
        Renumber residues incrementally to prevent problems with unusual residue numbering: 75A, 75B

        Parameters
        ----------
        prot: :class:moleculekit.moelcule.Molecule
            Molecule object of the protein to be renumbered

        """
        m = prot.renumberResidues(returnMapping=True)
        outpath, file = os.path.split(self.structure)
        file_name, file_extension = os.path.splitext(file)
        prot.write(os.path.join(outpath, f'{file_name}_renumbered{file_extension}'))
        m.to_csv(os.path.join(outpath, 'resid_map.csv'), index=False)

        logger.error(f'Your PDB structure has unusual residue numbering, please use the renumbered structure: {os.path.join(outpath, f"{file_name}_renumbered{file_extension}")} and the numbering provided in: {os.path.join(outpath, "resid_map.csv")} instead.')
        raise RuntimeError(f'Your PDB structure has unusual residue numbering, please use the renumbered structure: {os.path.join(outpath, f"{file_name}_renumbered{file_extension}")} and the numbering provided in: {os.path.join(outpath, "resid_map.csv")} instead.')

    def get_ff_params(self, lig_mol2: List = [], lig_frcmod: List = [], lig_rtf: List = [], lig_prm: List = []) -> \
            Union[Union[amber.defaultTopo, amber.defaultParam], Union[charmm.defaultTopo, charmm.defaultParam]]:
        """
        Returns the default force fields of HTMD for Amber and CHARMM:
        - Amber ff14sb
        - CHARMM36

        Parameters
        ----------
        lig_mol2: list
            List of .mol2 files for ligand to add (amber)
        lig_frcmod: list
            List of .frcmod parameter modification files for ligand to add (amber)
        lig_rtf: list
            List of .rtf topology files for ligand to add (charmm)
        lig_prm: list
            List of .prm parameter files for ligand to add (charmm)
        Returns
        -------
        List of parameter files that are necessary for building with HTMD.

        """
        if self.forcefield.startswith('amber'):
            return amber.defaultTopo() + lig_mol2, amber.defaultParam() + lig_frcmod
        elif self.forcefield.startswith('charmm'):
            return charmm.defaultTopo() + lig_rtf, charmm.defaultParam() + lig_prm
        else:
            logger.error(f'{self.forcefield} not implemented.')
            raise ValueError(f'{self.forcefield} not implemented.')

    def build_complex(self, sampling_pocket: str = 'GLY') -> NoReturn:
        """Builds proteins for force field computations with ligand inside.

        Force field computations need .param/.psf files for proteins or .frcmod/.rtf/.prm files for ligands, depending on the force field.
        These files describe a molecular system, as they contain all parameters.
        This procedure build different proteins that will be needed during the later design stages.
        Scaffold overview:
        - Glycin scaffold for sampling ligand poses
        - Native protein/ligand complex

        Writes the build systems into the '/my_project/scaffold/protein_params/' directory.

        Parameters
        ----------
        sampling_pocket: str
            What type of sampling pockets to build for rotamer and ligand pose sampling procedures
            In theory every pocket is possible, but only glycine and alanine are meanigful [default: GLY]
        """
        # Merge scaffold and ligand into one molecule
        prot = Molecule(self.structure)

        lig = Molecule(self.ligand_protonated)
        lig.set('segid', 'L')
        lig.set('chain', 'L')
        lig.set('resname', 'MOL')

        # mutate non-standard residues to prevent errors in building
        for non_standard, standard in special2standard.items():
            if non_standard in prot.resname:
                logger.info(f'Non-standard amino acid: {non_standard} found in protein, converting to: {standard}.')
                prot.mutateResidue(sel=f'resname {non_standard}', newres=standard)
        if 'UNK' in prot.resname:
            logger.error('Unknown amino acid in protein. No parametrization in system building possible. Remove it from protein or provide parameters.')
            raise RuntimeError('Unknown amino acid in protein. No parametrization in system building possible. Remove it from protein or provide parameters.')

        mol = Molecule(name='complex')
        mol.append(prot)
        mol.append(lig)

        # Build a native protein/ligand complex for minimization
        if not self.mutations:
            # Path to write prepared complex
            native_complex_path = os.path.join(self.work_dir, 'scaffold', self.forcefield, 'protein_params', 'native_complex')
            if not os.path.isfile(os.path.join(native_complex_path, 'structure.pdb')):
                logger.info('Build native complex.')
                self.build_ff(
                    structure=mol,
                    outdir=native_complex_path)
            else:
                logger.info('Native complex already built.')

        # Build all mutated scaffolds needed for energy calculations
        else:
            # Build glycine scaffold for ligand pose sampling
            self.ligand_sampling_pocket(struc=mol, sampling_pocket=sampling_pocket)

            # Build scaffolds with only a single mutation
            self.build_single_mut(struc=mol, sampling_pocket=sampling_pocket)

            # Build scaffolds containing pairwise mutations
            self.build_pair_mut(struc=mol)

    def build_single_mut(self, struc: Molecule, sampling_pocket: str) -> NoReturn:
        """
        Builds scaffolds with single mutations, where all other positions are mutated to either glycine or alanine

        Parameters
        ----------
        struc: :class: moleculekit.molecule.Molecule
            Object of molecule which will be mutated
        sampling_pocket: str
            What type of sampling pockets to build for rotamer and ligand pose sampling procedures
            In theory every pocket is possible, but only glycine and alanine are meaningful

        """
        # Build always again, since mutations could have changed
        for position in self.mutations:
            chain, resid = position['chain'], position['resid']
            for resname in position["mutations"]:
                outpath = os.path.join(self.work_dir, 'scaffold', self.forcefield, 'protein_params', f'{chain}_{resid}_{resname}')
                if not os.path.isfile(os.path.join(outpath, 'structure.pdb')):
                    logger.info(f'Build mutation: {chain}_{resid}_{resname}.')
                    sample_structure = struc.copy()
                    for other_position in self.mutations:
                        # Other position: mutate to glycine or alanine
                        if not other_position['chain'] == chain or not other_position['resid'] == resid:
                            sample_structure.mutateResidue(sel=f'chain {other_position["chain"]} and resid {other_position["resid"]}', newres=sampling_pocket)
                        # Same position
                        else:
                            # Residue is mutated otherwise Keep unmutated residues from the minimized structure
                            if not resname == sample_structure.get('resname', f'chain {chain} and resid {resid}')[0]:
                                sample_structure.mutateResidue(sel=f'chain {other_position["chain"]} and resid {other_position["resid"]}', newres=resname)
                    self.build_ff(
                        structure=sample_structure,
                        outdir=outpath
                    )
                else:
                    logger.info(f'Mutation {chain}_{resid}_{resname}, already built.')

    def ligand_sampling_pocket(self, struc: Molecule, sampling_pocket: str) -> NoReturn:
        """
        Build a glycine pocket, where all design positions are mutated to glycine to sample ligand poses

        Parameters
        ----------
        struc: :class: moleculekit.molecule.Molecule
            Object of molecule which will be mutated
        sampling_pocket: str
            What type of sampling pockets to build for rotamer and ligand pose sampling procedures
            In theory every pocket is possible, but only glycine and alanine are meanigful
        """
        # Build always again, since mutations could have changed
        sample_pocket_path = os.path.join(self.work_dir, 'scaffold', self.forcefield, 'protein_params', 'ligand_sampling_pocket')
        if not os.path.isfile(os.path.join(sample_pocket_path, 'structure.pdb')):
            logger.info('Build ligand sampling pocket.')
            sample_structure = struc.copy()
            for position in self.mutations:
                sample_structure.mutateResidue(sel=f'chain {position["chain"]} and resid {position["resid"]}', newres=sampling_pocket)

            self.build_ff(
                structure=sample_structure,
                outdir=sample_pocket_path
            )
        else:
            logger.info('Ligand sampling pocket already built.')

    def build_pair_mut(self, struc: Molecule) -> NoReturn:
        """ Builds a scaffold with two mutations inside.

        Parameters
        ----------
        struc: :class: moleculekit.molecule.Molecule
            Object of molecule which will be mutated
        """

        sorted_pair_permutations = create_pairs(self.mutations)

        for mutation_a, mutation_b in sorted_pair_permutations:
            # Get name of each residue from each pair
            chain_a, resid_a = mutation_a.split('_')[0:2]
            chain_b, resid_b = mutation_b.split('_')[0:2]
            resname_a = mutation_a.split('_')[2]
            resname_b = mutation_b.split('_')[2]
            mut_structure = struc.copy()
            scaffold_name = f'{chain_a}_{resid_a}_{resname_a}_{chain_b}_{resid_b}_{resname_b}'

            outpath = os.path.join(self.work_dir, 'scaffold', self.forcefield, 'protein_params', scaffold_name)

            if not os.path.isfile(os.path.join(outpath, 'structure.pdb')):
                logger.info(f'Build mutation combination: {scaffold_name}.')
                mut_structure.mutateResidue(sel=f'chain {chain_a} and resid {resid_a}', newres=resname_a)
                mut_structure.mutateResidue(sel=f'chain {chain_b} and resid {resid_b}', newres=resname_b)

                self.build_ff(
                    structure=mut_structure,
                    outdir=outpath)
            else:
                logger.info(f'Mutations {scaffold_name} already build.')

    def build_ff(self, structure: Molecule, outdir: str, renumber: bool = True) -> NoReturn:
        """ Basic force field building function.

        Parameters
        ----------
        structure: :class: moleculekit.molecule.Molecule
            Object of Molecule the system will be built for
        outdir: str
            Output path.
        renumber: bool
            Amber/CHARMM building renumbers resids consecutively instead of keeping the original resids. This option tries
            to renumber the output file according to the input structure. If it fails, the structure is renumbered from 0
            to X instead and a conversion table is written.
        """

        # Create the output directories
        os.makedirs(outdir, exist_ok=True)
        # Load default force field parameters
        if self.built_ligand_params:
            topo, param = self.get_ff_params(lig_mol2=self.built_ligand_params['mol2'],
                                             lig_frcmod=self.built_ligand_params['frcmod'],
                                             lig_rtf=self.built_ligand_params['rtf'],
                                             lig_prm=self.built_ligand_params['prm'])
        else:
            topo, param = self.get_ff_params()

        segids = np.unique(structure.segid)
        caps = {}

        if self.forcefield.startswith('amber'):
            for i in segids:
                if i == 'L':
                    continue
                caps[i] = ['none', 'none']

            # Only create the input files with htmd, since we need to modify the tleap.in file
            amber.build(
                mol=structure,
                ff=['leaprc.gaff2', 'leaprc.protein.ff14SB'],
                topo=topo,
                param=param,
                ionize=False,
                outdir=outdir,
                caps=caps,
                execute=False
            )
            if os.path.isfile(os.path.join(outdir, 'tleap.in')):
                with open(os.path.join(outdir, 'tleap.in'), 'r') as input_file:
                    with open(os.path.join(outdir, 'tleap_corrected.in'), 'w') as output_file:
                        for line in input_file:
                            # We need to remove this line since it moves the protein
                            if not line.startswith('setBox mol "vdw"'):
                                output_file.write(line)
                os.rename(os.path.join(outdir, 'tleap_corrected.in'), os.path.join(outdir, 'tleap.in'))

                built_command = [self.tleap, '-f', './tleap.in']

                logger.info('Starting the build.')
                os.chdir(outdir)
                subprocess.run(built_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                logger.info('Finished building.')

                if os.path.exists(os.path.join(outdir, "structure.crd")) and os.path.exists(os.path.join(outdir, 'structure.prmtop')):
                    try:
                        molbuilt = Molecule(os.path.join(outdir, 'structure.prmtop'))
                        molbuilt.read(os.path.join(outdir, 'structure.crd'))
                    except Exception as e:
                        raise RuntimeError(
                            f'Failed at reading structure.prmtop/structure.crd due to error: {e}'
                        )
                else:
                    raise RuntimeError(
                        f'No structure pdb/prmtop file was generated. Check {os.path.join(outdir, "leap.log")} for errors in building.')
            else:
                logger.error('Failed to create tleap input file.')
                raise FileNotFoundError('Failed to create tleap input file.')

            molbuilt.write(os.path.join(outdir, 'structure.pdb'))
            bmol = Molecule(os.path.join(outdir, 'structure.pdb'))

        elif self.forcefield.startswith('charmm'):
            for i in segids:
                if i == 'L':
                    continue
                caps[i] = ['first none', 'last none']
            bmol = charmm.build(
                mol=structure,
                topo=topo,
                param=param,
                ionize=False,
                outdir=outdir,
                caps=caps,
                psfgen=path(bin_file='psfgen')
            )
            _ = fix_parameters(os.path.join(outdir, 'parameters'))
        else:
            raise ValueError('No valid force field defined!')

        bmol = autoSegment2(bmol, _logger=False)

        bmol.set('segid', 'L', f'resname MOL')
        bmol.set('chain', 'L', f'resname MOL')

        segids = list(filter(None, np.unique(bmol.segid)))
        for segid in segids:
            chain = np.unique(structure.get('chain', sel=f'segid {segid}'))[0]
            bmol.set('chain', chain, f'segid {segid}')

        if renumber:
            try:
                bmol.renumberResidues()
                renum_resids = [k for k, g in itertools.groupby(bmol.resid)]
                old_resids = [k for k, g in itertools.groupby(structure.resid)]
                resmap = dict(zip(renum_resids, old_resids))
                count_dups = [sum(1 for _ in group) for _, group in itertools.groupby(bmol.resid)]
                assert len(old_resids) == len(renum_resids)
                assert np.unique(structure.chain).shape[0] == np.unique(bmol.chain).shape[0]
                count = 0
                for i, e in enumerate(count_dups):
                    bmol.resid[count:count + e] = resmap[i]
                    count += e
            except AssertionError:
                logger.warning(f'Check the resid_map.csv under {outdir} to get an overview of the new resids.')
                m = structure.renumberResidues(returnMapping=True)
                m.to_csv(os.path.join(outdir, 'resid_map.csv'), index=False)
        bmol.write(os.path.join(outdir, 'structure.pdb'))

    def parameterize_ligand(self, addHs: bool = True) -> NoReturn:
        """ Ligand preparation procedure

        Steps:
        1. Add hydrogens with OpenBabel
        2. Assign force field parameters and partial charges for GAFF2 with ANTECHAMBER
        3. Generate .frcmod with PARMCHK2

        Or:

        2. Generate .prm/.rtf with MATCH-TYPER

        Parameters
        ---------

        addHs: bool
            Whether to add hydrogen atoms [default: True]
        """

        if not os.path.isfile(self.structure):
            logger.error(f'Could not open: {self.structure}.')
            raise FileNotFoundError(f'Could not open: {self.structure}.')

        ligand_mol2 = os.path.join(self.built_ligand_params['params_folder'], 'ligand.mol2')
        ligand_pdb = os.path.join(self.built_ligand_params['params_folder'], 'ligand.pdb')
        os.makedirs(self.built_ligand_params['params_folder'], exist_ok=True)
        os.chdir(self.built_ligand_params['params_folder'])

        # Remove all protons and convert to mol2
        ligand = Molecule(self.structure)
        if addHs:
            ligand.remove('element H', _logger=False)
        ligand.write(ligand_mol2)

        if addHs:
            # Protonate with Obabel
            logger.info(f'Adding hydrogen atoms to the ligand according to pH: {self.ph}.')
            protonate_command = [self.obabel, '-i', 'mol2', ligand_mol2, '-o', 'mol2', '-O', ligand_mol2, '--p', str(self.ph)]
            process = subprocess.Popen(protonate_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
            if process.returncode != 0:
                logger.error(f'Obabel failed with the following exception: {stderr.decode("ascii")}')
                raise RuntimeError(f'Obabel failed with the following error: {stderr.decode("ascii")}')

        # Calculate netcharge for parameterization
        ligand = Molecule(ligand_mol2)
        netcharge = str(round(np.sum(ligand.charge)))
        ligand.write(ligand_pdb)

        # Change the residue name to MOL, assign bond information and partial charges
        conv_command = [self.antechamber, '-at', 'sybyl', '-c', 'gas', '-rn', 'MOL', '-nc', netcharge,
                        '-i', ligand_pdb, '-o', ligand_mol2, '-fi', 'pdb', '-fo', 'mol2', '-j', '4']

        process = subprocess.Popen(conv_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            logger.error(f'Antechamber failed with the following exception: {stderr.decode("ascii")}')
            raise RuntimeError(f'Antechamber failed with the following error: {stderr.decode("ascii")}')

        # Copy prepared ligand
        shutil.copy(ligand_mol2, os.path.abspath(os.path.join('..', 'ligand.mol2')))

        # Convert atomtypes to GAFF2 for AMBER
        if self.forcefield.startswith('amber'):
            ff = 'GAFF2'
            logger.info(f'Parameterize ligand for {ff}.')
            antechamber_command = [
                self.antechamber,
                '-at',
                'gaff2',
                '-c',
                'gas',
                '-rn',
                'MOL',
                '-nc',
                netcharge,
                '-fi',
                'pdb',
                '-i',
                ligand_pdb,
                '-fo',
                'mol2',
                '-o',
                ligand_mol2,
                '-j',
                '4']

            process = subprocess.Popen(antechamber_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
            if process.returncode != 0:
                logger.error(f'Antechamber failed with the following error: {stderr.decode("ascii")}')
                raise RuntimeError(f'Antechamber failed with the following error: {stderr.decode("ascii")}')

            # Create .frcmod file
            parmchk2_command = [
                self.parmchk2,
                '-f',
                'mol2',
                '-s',
                'gaff2',
                '-i',
                ligand_mol2,
                '-o',
                'ligand.frcmod']
            process = subprocess.Popen(parmchk2_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
            if process.returncode != 0:
                logger.error(f'Parmchk2 failed with the following error: {stderr.decode("ascii")}')
                raise RuntimeError(f'Parmchk2 failed with the following error: {stderr.decode("ascii")}')

        elif self.forcefield.startswith('charmm'):
            ff = 'CGenFF_2b6'
            logger.info(f'Parameterize ligand for {ff}.')
            match_cmd = [
                self.match,
                '-forcefield',
                'top_all36_cgenff_new',
                ligand_mol2]

            process = subprocess.Popen(match_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
            if process.returncode != 0:
                logger.error(f'Match-typer failed with the following error: {stderr.decode("ascii")}')
                raise RuntimeError(f'Match-typer failed with the following error: {stderr.decode("ascii")}')

        logger.info('Ligand parametrization was successful.')
        os.chdir(self.work_dir)

    def minimize_structure(self, structure_path: str, cuda: bool = False, restraint_bb: bool = True) -> NoReturn:
        """
        Minimization method utilizing OpenMM.
        The structures will be initialized with the corresponding force
        field and then be minimized until the energy converges.

        Parameters
        ----------
        structure_path: str
            Path to built complex
        cuda: bool
            Minimization on GPU. [default: False]
        restraint_bb: bool
            Applies a restraint on the backbone. [default: True]
        """

        if self.forcefield.startswith('amber'):
            # Load Amber input files
            structure, params = load_ff_parameters(structure_path=structure_path, forcefield=self.forcefield)
            structure_prm = app.AmberPrmtopFile(os.path.join(structure_path, 'structure.prmtop'))
            inpcrd = app.AmberInpcrdFile(os.path.join(structure_path, 'structure.crd'))
            system = structure_prm.createSystem()
        elif self.forcefield.startswith('charmm'):
            # Load Charmm input files
            structure, params = load_ff_parameters(structure_path=structure_path, forcefield=self.forcefield)
            structure_prm = parmed.charmm.CharmmPsfFile(os.path.join(structure_path, 'structure.psf'))
            inpcrd = app.PDBFile(os.path.join(structure_path, 'structure.pdb'))
            system = structure_prm.createSystem(params)
        else:
            logger.error('Force field not supported.')
            raise ValueError('Force field not supported.')
        metric_rmsd = MetricRmsd(structure.copy(), 'all', pbc=False)

        # Standard values for the integrator and tolerance constraint
        temperature = self.temperature * unit.kelvin
        friction = 1.0 / unit.picoseconds
        error_tolerance = 2.0 * unit.femtoseconds
        if not cuda:
            # Default value
            tolerance = 10.0 * unit.kilojoule/unit.mole
        else:
            # High tolerance so the CPU only pre-minimizes
            tolerance = 1e6

        # Prepare System and Integrator
        cpu_integrator = mm.VariableLangevinIntegrator(
            temperature,
            friction,
            error_tolerance
        )

        if restraint_bb:
            # Applies an external force on backbone atoms
            # This allows the backbone to stay rigid, while severe clashes can still be resolved
            force = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
            force.addGlobalParameter("k", 5.0 * unit.kilocalories_per_mole / unit.angstroms ** 2)
            force.addPerParticleParameter("x0")
            force.addPerParticleParameter("y0")
            force.addPerParticleParameter("z0")
            for i, atom_crd in enumerate(inpcrd.positions):
                if structure.name[i] in ('N', 'CA', 'C', 'O'):
                    force.addParticle(i, atom_crd.value_in_unit(unit.nanometers))
            system.addForce(force)

        # Setup CPU minimization
        platform = mm.Platform.getPlatformByName('CPU')
        properties = None
        cpu_min = app.Simulation(structure_prm.topology, system, cpu_integrator, platform, properties)
        cpu_min.context.setPositions(inpcrd.positions)
        pre_min_state = cpu_min.context.getState(getEnergy=True, getForces=True)
        pre_min_state_nrg = pre_min_state.getPotentialEnergy().in_units_of(unit.kilocalories_per_mole)
        logger.info(f'Starting potential energy: {pre_min_state_nrg}.')
        # The energy will only be minimized to the set tolerance value
        cpu_min.minimizeEnergy(tolerance=tolerance)
        position = cpu_min.context.getState(getPositions=True).getPositions()
        post_min_state = cpu_min.context.getState(getEnergy=True, getForces=True)

        if cuda:
            # Get the pre-minimized state
            min_coords = cpu_min.context.getState(getPositions=True)
            # Setup GPU minimizer and integrator
            platform = mm.Platform.getPlatformByName('CUDA')
            properties = {'CudaPrecision': 'mixed'}
            gpu_integrator = mm.VariableLangevinIntegrator(
                temperature,
                friction,
                error_tolerance
            )
            gpu_min = app.Simulation(structure_prm.topology, system, gpu_integrator, platform, properties)
            gpu_min.context.setPositions(min_coords.getPositions())
            # Minimize until convergence
            gpu_min.minimizeEnergy()
            position = gpu_min.context.getState(getPositions=True).getPositions()
            post_min_state = gpu_min.context.getState(getEnergy=True, getForces=True)
        post_min_state_nrg = post_min_state.getPotentialEnergy().in_units_of(unit.kilocalories_per_mole)
        logger.info(f'Ending potential energy: {post_min_state_nrg}.')

        app.PDBFile.writeFile(structure_prm.topology, position, open(self.built_scaffold, 'w'), keepIds=True)

        # The OpenMM output PDB has usally messed up chains, segids etc.
        # So we just set the minimized coordinates to the input structure
        prot_min = Molecule(self.built_scaffold)
        rmsd = float(metric_rmsd.project(prot_min))
        logger.info(f'RMSD between minimized and unminimized structure: {str(round(rmsd, 4))} Å.')
        structure.coords = prot_min.coords
        structure.write(self.built_scaffold, sel='not segid L')
