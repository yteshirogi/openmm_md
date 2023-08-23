import openmm as mm
import openmm.app as app
from openmm import unit


def generate_system(complex_topology, forcefield, complex_positions):
    """
    Generate an OpenMM System object from a topology, forcefield and positions.

    Parameters
    ----------
    complex_topology: openmm.app.Topology
        Topology of the system.
    forcefield: openmm.app.Forcefield
        Forcefield of the system.
    complex_positions: openmm.unit.Quantity
        Positions of the system.

    Returns
    -------
    simulation: openmm.app.Simulation
        Simulation object.
    """
    modeller = app.Modeller(complex_topology, complex_positions)
    modeller.addSolvent(
        forcefield,
        padding=1.0*unit.nanometers,
        ionicStrength=0.15*unit.molar
    )
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0*unit.nanometers,
        constraints = app.HBonds,
        rigidWater=True,
        ewaldErrorTolerance=0.0005
    )
    system.addForce(
        mm.MonteCarloBarostat(1.0*unit.atmospheres, 300*unit.kelvin, 25)
    )
    integrator = mm.LangevinIntegrator(
        300*unit.kelvin,
        1.0/unit.picoseconds,
        2.0*unit.femtoseconds
    )
    simulation = app.Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)

    return simulation
