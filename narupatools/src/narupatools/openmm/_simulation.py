# This file is part of narupatools (https://github.com/alexjbinnie/narupatools).
# Copyright (c) Alex Jamieson-Binnie. All rights reserved.
#
# narupatools is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# narupatools is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with narupatools.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import annotations

import contextlib
from typing import Generator, Iterator, Mapping, MutableMapping, Optional

from infinite_sets import InfiniteSet
from narupa.trajectory import FrameData
from openmm import (
    Context,
    Integrator,
    LocalEnergyMinimizer,
    OpenMMException,
    Platform,
    System,
)
from openmm.app import Simulation, Topology

from narupatools.frame import FrameSource
from narupatools.override import override

from ._converter import openmm_simulation_to_frame
from ._serializer import deserialize_simulation


class OpenMMSimulation(FrameSource):
    """
    Alternative representation of an OpenMM simulation.

    An OpenMM simulation as defined by OpenMM is a collection of a system, topology, integrator
    and context. In OpenMM, these objects are represented as C++ objects, but the simulation itself
    is a pure Python object. However, it is not easy to manipulate the parts of the simulation without
    having to construct a whole new OpenMM simulation object.

    This alternative implementation acts similarly, but allows the system, integrator and topology to
    be easily modified. It also provides some more pythonic ways to interact with your system.
    """

    @override(FrameSource.get_frame)
    def get_frame(  # noqa: D102
        self, *, fields: InfiniteSet[str], existing: Optional[FrameData] = None
    ) -> FrameData:

        frame = existing
        if not frame:
            frame = FrameData()
        openmm_simulation_to_frame(self, fields=fields, existing=frame)  # type: ignore
        return frame

    def __init__(
        self,
        *,
        system: System,
        topology: Topology,
        integrator: Integrator,
        context: Optional[Context] = None,
        platform: Optional[Platform] = None,
    ):
        self._system = system
        self._topology = topology
        self._integrator = integrator
        if context:
            self._context = context
        else:
            self._reset_context()
        if platform:
            self._platform = platform
        else:
            self._platform = self._context.getPlatform()

    def _reset_context(self, recreate: bool = False) -> None:
        if recreate:
            state = self._context.getState(getPositions=True, getVelocities=True)
            self._context = Context(self.system, self.integrator)
            self._context.setPositions(state.getPositions())  # type: ignore
            self._context.setVelocities(state.getVelocities())  # type: ignore
        else:
            self._context.reinitialize(preserveState=True)

    @contextlib.contextmanager
    def modify(self) -> Generator[OpenMMSimulation, None, None]:
        """
        Context manager for ensuring the context is reset after modifications are made.

        Use this when you intend to modify the system, integrator or topology 'in-place'::

           with simulation.modify():
               simulation.system.addForce(...)
               ...

           # the context will automatically be regenerated when the block is finished.
        """
        yield self
        self._reset_context()

    @property
    def system(self) -> System:
        """System describing the forces and masses of the simulation."""
        return self._system

    @system.setter
    def system(self, value: System) -> None:
        self._system = value
        self._reset_context()

    @property
    def integrator(self) -> Integrator:
        """Integrator used to propagate the simulation forward."""
        return self._integrator

    @integrator.setter
    def integrator(self, value: Integrator) -> None:
        self._integrator = value
        self._reset_context(True)

    @property
    def topology(self) -> Topology:
        """Topology describing the bonds and residues of the simulation."""
        return self._topology

    @topology.setter
    def topology(self, value: Topology) -> None:
        self._topology = value
        self._reset_context()

    @property
    def context(self) -> Context:
        """Context referring to current calculations and state of the simulation."""
        return self._context

    @classmethod
    def from_simulation(cls, simulation: Simulation) -> OpenMMSimulation:
        """
        Create a simulation that is identical to the given OpenMM simulation.

        This does not copy the system, integrator, topology or context, but merely creates
        an OpenMMSimulation that refers to the same object.
        """
        return cls(
            system=simulation.system,
            topology=simulation.topology,
            integrator=simulation.integrator,
            context=simulation.context,
        )

    @classmethod
    def from_xml_file(cls, filename: str) -> OpenMMSimulation:
        """Create a simulation from an XML filename."""
        with open(filename) as file:
            return cls.from_simulation(deserialize_simulation(file.read()))

    def minimize(
        self, tolerance: float = 1, max_iterations: Optional[int] = None
    ) -> None:
        """
        Minimize the energy using a L-BFGS minimizer.

        :param tolerance: Desired energy convergence criterea in kilojoules per mole.
        :param max_iterations: Maximum iterations to run, or None to run as many as required.
        """
        if not max_iterations:
            max_iterations = 0
        LocalEnergyMinimizer.minimize(self._context, tolerance, max_iterations)

    def run(self, steps: int) -> None:
        """Run the simulation for the given number of steps."""
        self.integrator.step(steps)

    @property
    def global_parameters(self) -> OpenMMParametersView:
        """View of the global parameters of the system."""
        return OpenMMParametersView(self)

    def create_checkpoint(self) -> str:
        """Generate a checkpoint representing a snapshot of the system."""
        return self._context.createCheckpoint()

    def load_checkpoint(self, checkpoint: str) -> None:
        """Reinitialize the simulation from a checkpoint."""
        self._context.loadCheckpoint(checkpoint)
        self._context.reinitialize(preserveState=True)


class OpenMMParametersView(MutableMapping[str, float]):
    """
    View of the global parameters of an OpenMM simulation.

    This is a dynamic view which updates when the underlying simulation is changed. It acts as a
    mutable dictionary, except keys may not be added or deleted.
    """

    def __init__(self, simulation: OpenMMSimulation):
        self._simulation = simulation

    def __setitem__(self, k: str, v: float) -> None:
        try:
            return self._simulation._context.setParameter(k, v)
        except OpenMMException as e:
            raise KeyError(f"No parameter with name {k}") from e

    def __delitem__(self, v: str) -> None:
        raise NotImplementedError("Cannot delete global parameter.")

    def __getitem__(self, k: str) -> float:
        try:
            return self._simulation._context.getParameter(k)
        except OpenMMException as e:
            raise KeyError(f"No parameter with name {k}") from e

    def __len__(self) -> int:
        return len(self.__parameters)

    def __iter__(self) -> Iterator[str]:
        yield from self.__parameters.keys()

    @property
    def __parameters(self) -> Mapping[str, float]:
        return self._simulation._context.getParameters()

    def __str__(self) -> str:
        return f"<OpenMMParametersView {len(self)} parameter(s)>"
