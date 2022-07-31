import os
import logging

from openmm import app
import openmm as mm
from simtk import unit
import parmed
from typing import NoReturn
from moleculekit.molecule import Molecule
from moleculekit.projections.metricrmsd import MetricRmsd

from pocketoptimizer.utility.utils import load_ff_parameters
from pocketoptimizer.utility.molecule_types import _BB_ATOMS

logger = logging.getLogger(__name__)


def minimize_structure(structure_path: str, forcefield: str, output_pdb: str, cuda: bool = False, restraint_bb: bool = True, temperature: float = 300.0) -> NoReturn:
    """
    Minimization method utilizing OpenMM.
    The structures will be initialized with the corresponding force
    field and then be minimized until the energy converges.

    Parameters
    ----------
    structure_path: str
        Path to the built structures.
    forcefield: str
        Chosen force field (i.e. charmm36).
    output_pdb:
        Output name of the minimized structure.
    cuda: bool
        Minimization on GPU. [default: False]
    restraint_bb: bool
        Applies a restraint on the backbone. [default: True]
    temperature: float
        temperature of the system in Kelvin [default: 300.0]
    """

    if forcefield.startswith('amber'):
        # Load Amber input files
        structure, params = load_ff_parameters(structure_path=structure_path, forcefield=forcefield)
        structure_prm = app.AmberPrmtopFile(os.path.join(structure_path, 'structure.prmtop'))
        inpcrd = app.AmberInpcrdFile(os.path.join(structure_path, 'structure.crd'))
        system = structure_prm.createSystem()
    elif forcefield.startswith('charmm'):
        # Load Charmm input files
        structure, params = load_ff_parameters(structure_path=structure_path, forcefield=forcefield)
        structure_prm = parmed.charmm.CharmmPsfFile(os.path.join(structure_path, 'structure.psf'))
        inpcrd = app.PDBFile(os.path.join(structure_path, 'structure.pdb'))
        system = structure_prm.createSystem(params)
    else:
        logger.error('Force field not supported.')
        raise ValueError('Force field not supported.')
    metric_rmsd = MetricRmsd(structure.copy(), 'all', pbc=False)

    # Standard values for the integrator and tolerance constraint
    temperature = temperature * unit.kelvin
    friction = 1.0 / unit.picoseconds
    error_tolerance = 2.0 * unit.femtoseconds
    distance_tolerance = 0.00001
    if not cuda:
        # Default value
        tolerance = 10 * unit.kilojoule/unit.mole
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
            if structure.name[i] in _BB_ATOMS and not structure.segid[i] == 'L':
                force.addParticle(i, atom_crd.value_in_unit(unit.nanometers))
        system.addForce(force)

    # Setup CPU minimization
    cpu_integrator.setConstraintTolerance(distance_tolerance)
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
        gpu_integrator.setConstraintTolerance(distance_tolerance)
        gpu_min = app.Simulation(structure_prm.topology, system, gpu_integrator, platform, properties)
        gpu_min.context.setPositions(min_coords.getPositions())
        # Minimize until convergence
        gpu_min.minimizeEnergy()
        position = gpu_min.context.getState(getPositions=True).getPositions()
        post_min_state = gpu_min.context.getState(getEnergy=True, getForces=True)
    post_min_state_nrg = post_min_state.getPotentialEnergy().in_units_of(unit.kilocalories_per_mole)
    logger.info(f'Ending potential energy: {post_min_state_nrg}.')

    app.PDBFile.writeFile(structure_prm.topology, position, open(output_pdb, 'w'), keepIds=True)

    # The OpenMM output PDB has usally messed up chains, segids etc.
    # So we just set the minimized coordinates to the input structure
    prot_min = Molecule(output_pdb)
    rmsd = float(metric_rmsd.project(prot_min))
    logger.info(f'RMSD between minimized and unminimized structure: {str(round(rmsd, 4))} Å.')
    structure.coords = prot_min.coords
    structure.remove('segid L', _logger=False)
    structure.write(output_pdb)
