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

from typing import Any, Dict, Optional, Protocol

import numpy as np
from ase.atoms import Atoms

from narupatools.ase._units import UnitsASE
from narupatools.core import UnitsNarupa
from narupatools.imd import Interaction
from narupatools.physics.typing import Vector3Array

from ._constraint import ASEEnergyConstraint, ASEMomentaConstraint, ASETorqueConstraint
from ._null_constraint import NullConstraint

_NarupaToASE = UnitsNarupa >> UnitsASE
_ASEToNarupa = UnitsASE >> UnitsNarupa


class ASEAtomsWrapper(Protocol):
    """Wrappar around an ASE atoms object."""

    @property
    def atoms(self) -> Atoms:
        """Wrapped ASE Atoms."""
        ...


class _ASEAtomsWrapper(ASEAtomsWrapper):
    def __init__(self, atoms: Atoms):
        self._atoms = atoms

    @property
    def atoms(self) -> Atoms:
        return self._atoms


class InteractionConstraint(ASEEnergyConstraint, ASEMomentaConstraint, ASETorqueConstraint):
    """
    An ASE constraint that applies an iMD force.

    This conflates both an implementation of a narupatools interaction and an ASE
    constraint.
    """

    def __init__(
        self,
        *,
        interaction: Interaction,
    ):
        """
        Create an ASE constraint that will apply an interactive force.

        :param dynamics: Dynamics this interaction is attached to.
        :param interaction: Initial parameters of the interaction.
        """
        self.interaction = interaction

    def adjust_positions(  # noqa: D102
        self, atoms: Atoms, positions: np.ndarray, /
    ) -> None:
        # Assume all interactions depend on positions
        self.interaction.calculate_forces_and_energy()

    def adjust_momenta(  # noqa: D102
        self, atoms: Atoms, momenta: np.ndarray, /
    ) -> None:
        # Assume no interactions depend on velocities
        # When they do, this should conditionally invalidate the cache
        self.interaction.calculate_forces_and_energy()

    def adjust_forces(self, atoms: Atoms, forces: np.ndarray, /) -> None:  # noqa: D102
        forces[self.interaction.particle_indices] += self.interaction.forces * _NarupaToASE.force

    def adjust_potential_energy(self, /, atoms: Atoms) -> float:  # noqa: D102
        return self.interaction.potential_energy * _NarupaToASE.energy

    def adjust_torques(self, atoms: Atoms, torques: np.ndarray, /) -> None:
        torques[self.interaction.particle_indices] += self.interaction.torques * _NarupaToASE.force

    def __deepcopy__(self, memodict: Dict[Any, Any] = None) -> Any:
        # Deep copying an IMD constraint destroys it. This is because the interaction itself
        # has a reference to the atoms object, and this leads to deep copying of the atoms
        # object and its calculator, which may not be designed to be deep copied.
        return NullConstraint()
