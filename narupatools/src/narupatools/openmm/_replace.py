import copy

from openmm import Integrator
from openmm.app import Simulation


def replace_integrator(simulation: Simulation, integrator: Integrator) -> Simulation:
    new_simulation = Simulation(simulation.topology, simulation.system, integrator, platform=simulation.context.getPlatform())

    state = simulation.context.getState(getPositions=True, getVelocities=True)
    new_simulation.context.setPositions(state.getPositions(asNumpy=True))
    new_simulation.context.setVelocities(state.getVelocities(asNumpy=True))

    return new_simulation


def recreate_simulation(simulation: Simulation):
    new_simulation = Simulation(copy.deepcopy(simulation.topology), copy.deepcopy(simulation.system), copy.deepcopy(simulation.integrator), platform=simulation.context.getPlatform())

    state = simulation.context.getState(getPositions=True, getVelocities=True)
    new_simulation.context.setPositions(state.getPositions(asNumpy=True))
    new_simulation.context.setVelocities(state.getVelocities(asNumpy=True))

    return new_simulation
