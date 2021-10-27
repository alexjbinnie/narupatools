from contextlib import contextmanager

from openmm import System, Integrator, Context, LocalEnergyMinimizer
from openmm.app import Topology, Simulation


class OpenMMSimulation:

    def __init__(self, *, system: System, topology: Topology, integrator: Integrator, context: Context):
        self._system = system
        self._topology = topology
        self._integrator = integrator
        self._context = context

    @classmethod
    def from_simulation(cls, simulation: Simulation):
        return cls(system=simulation.system, topology=simulation.topology, integrator=simulation.integrator,
                   context=simulation.context)

    def _reset_context(self, recreate=False):
        if recreate:
            state = self._context.getState(getPositions=True, getVelocities=True)
            self._context = Context(self.system, self.integrator, self._context.getPlatform())
            self._context.setPositions(state.getPositions())
            self._context.setVelocities(state.getVelocities())
        else:
            self._context.reinitialize(preserveState=True)

    @property
    def system(self):
        return self._system

    @system.setter
    def system(self, value):
        self._system = value
        self._reset_context()

    @property
    def topology(self):
        return self._topology

    @topology.setter
    def topology(self, value):
        self._topology = value
        self._reset_context()

    @property
    def integrator(self):
        return self._integrator

    @integrator.setter
    def integrator(self, value):
        self._integrator = value
        self._reset_context(recreate=True)

    @property
    def context(self):
        return self._context

    @contextmanager
    def modify(self):
        yield self
        self._reset_context()

    def run(self, steps):
        self._integrator.step(steps)

    def minimize(self, tolerance, max_iterations=None):
        if not max_iterations:
            max_iterations = 0
        LocalEnergyMinimizer.minimize(self._context, tolerance, max_iterations)

    def create_checkpoint(self):
        return self._context.createCheckpoint()

    def load_checkpoint(self, checkpoint):
        self._context.loadCheckpoint(checkpoint)
        self._context.reinitialize(preserveState=False)
