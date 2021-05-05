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

"""An ASE constraint that triggers events when things change."""

import numpy as np
from ase.atoms import Atoms
from narupa.utilities.event import Event

from narupatools.core.event import EventListener
from .constraint import ASECellConstraint, ASEMomentaConstraint


class ASEObserver(ASECellConstraint, ASEMomentaConstraint):
    """
    ASE Constraint which listens for a change in its Atoms object.

    It is not advisable to act immediately when one of the events are triggered, as
    there may be subsequent constraints that will modify the value.
    """

    def __init__(self) -> None:
        """Create a new `ASEObserver`."""
        self._on_set_positions = Event()
        self._on_set_cell = Event()
        self._on_set_momenta = Event()

    @property
    def on_set_positions(self) -> EventListener:
        """Event triggered when the system has its positions modified."""
        return self._on_set_positions

    @property
    def on_set_momenta(self) -> EventListener:
        """Event triggered when the system has its momenta (and velocities) modified."""
        return self._on_set_momenta

    @property
    def on_set_cell(self) -> EventListener:
        """Event triggered when the system has its cell vectors modified."""
        return self._on_set_cell

    def adjust_positions(
        self, atoms: Atoms, positions: np.ndarray, /
    ) -> None:  # noqa: D102
        self._on_set_positions.invoke()

    def adjust_forces(self, atoms: Atoms, forces: np.ndarray, /) -> None:  # noqa: D102
        pass

    def adjust_cell(self, atoms: Atoms, cell: np.ndarray, /) -> None:  # noqa: D102
        self._on_set_cell.invoke()

    def adjust_momenta(
        self, atoms: Atoms, momenta: np.ndarray, /
    ) -> None:  # noqa: D102
        self._on_set_momenta.invoke()
