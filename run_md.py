from pathlib import Path

import mdtraj as md
import openmm.app as app
from openmm import unit


def run_simulation(simulation, steps=500000, write_interval=1000, log_interval=1000, DATA=Path('./data')):
    """
    Run a simulation.

    Parameters
    ----------
    simulation: openmm.app.Simulation
        Simulation object.
    """
    print('Performing energy minimization...')
    simulation.minimizeEnergy()
    with open(DATA / "topology.pdb", "w") as pdb_file:
        app.PDBFile.writeFile(
            simulation.topology,
            simulation.context.getState(
                getPositions=True,
                enforcePeriodicBox=True
            ).getPositions(),
            file=pdb_file,
            keepIds=True,
        )

    print('Equilibrating...')
    simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
    simulation.step(50000)

    print('Running Production...')
    simulation.reporters.append(
        md.reporters.XTCReporter(
            file=str(DATA / "trajectory.xtc"),
            reportInterval=write_interval
        )
    )
    simulation.reporters.append(
        app.StateDataReporter(
            str(DATA / "md.log"),
            log_interval,
            step=True,
            potentialEnergy=True,
            temperature=True,
            progress=True,
            remainingTime=True,
            speed=True,
            totalSteps=steps,
            separator="\t",
        )
    )
    simulation.reporters.append(
        app.CheckpointReporter(
            str(DATA / 'checkpoint.chk'),
            write_interval
        )
    )

    simulation.currentStep = 0
    simulation.step(steps)
