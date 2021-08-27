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

from narupatools.app._session import Session
from narupatools.core.broadcastable import Broadcastable, Broadcaster
from narupatools.core.dynamics import SimulationDynamics
from narupatools.imd._feature import InteractionFeature


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
            self.imd.add_source(broadcaster.shared_state.interactions.snapshot)

    def end_broadcast(self, broadcaster: Broadcaster) -> None:  # noqa: D102
        if isinstance(broadcaster, Session):
            self.imd.remove_source(broadcaster.shared_state.interactions.snapshot)
