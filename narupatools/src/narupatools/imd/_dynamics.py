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

from narupatools.app._session import Session, Broadcastable
from narupatools.core.dynamics import SimulationDynamics
from narupatools.imd._feature import InteractionFeature
from narupatools.override import override


class InteractiveSimulationDynamics(
    SimulationDynamics, Broadcastable, metaclass=ABCMeta
):
    """Base class for any simulation dynamics which supports IMD."""

    @property
    @abstractmethod
    def imd(self) -> InteractionFeature:
        """Access to the interactions currently being applied to the simulation."""
        raise NotImplementedError

    @override
    def start_broadcast(self, session: Session) -> None:  # noqa: D102
        self.imd.add_source(session.shared_state.interactions.snapshot)

    @override
    def end_broadcast(self, session: Session) -> None:  # noqa: D102
        self.imd.remove_source(session.shared_state.interactions.snapshot)
