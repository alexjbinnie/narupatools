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

"""Base class for simulations that support IMD."""

from abc import ABCMeta, abstractmethod

from narupatools.app.session import Session
from narupatools.core.broadcastable import Broadcastable, Broadcaster
from narupatools.core.dynamics import SimulationDynamics
from narupatools.imd.feature import InteractionFeature
from narupatools.physics.typing import ScalarArray, Vector3Array


class InteractiveSimulationDynamics(
    SimulationDynamics, Broadcastable, metaclass=ABCMeta
):
    """Base class for any simulation dynamics which supports IMD."""

    @property
    @abstractmethod
    def imd(self) -> InteractionFeature:
        """Access to the interactions currently being applied to the simulation."""
        raise NotImplementedError()

    def start_broadcast(self, broadcaster: Broadcaster) -> None:  # noqa: D102
        if isinstance(broadcaster, Session):
            self.imd.add_dynamic_interactions_source(broadcaster.server.imd)

    def end_broadcast(self, broadcaster: Broadcaster) -> None:  # noqa: D102
        if isinstance(broadcaster, Session):
            self.imd.remove_dynamic_interactions_source(broadcaster.server.imd)

    @property
    @abstractmethod
    def positions(self) -> Vector3Array:
        """Positions of particles in nanometers."""
        raise NotImplementedError()

    @property
    @abstractmethod
    def velocities(self) -> Vector3Array:
        """Velocities of particles in nanometers per picosecond."""
        raise NotImplementedError()

    @property
    @abstractmethod
    def forces(self) -> Vector3Array:
        """Forces on particles in kilojoules per mole per nanometer."""
        raise NotImplementedError()

    @property
    @abstractmethod
    def masses(self) -> ScalarArray:
        """Masses of particles in daltons."""
        raise NotImplementedError()

    @property
    @abstractmethod
    def kinetic_energy(self) -> float:
        """Kinetic energy in kilojoules per mole."""
        raise NotImplementedError()

    @property
    @abstractmethod
    def potential_energy(self) -> float:
        """Potential energy in kilojoules per mole."""
        raise NotImplementedError()
