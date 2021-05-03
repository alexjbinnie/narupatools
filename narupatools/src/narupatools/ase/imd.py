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

"""Implementation of IMD for ASE as dynamic restraints."""

from __future__ import annotations

from typing import Protocol

from ase import Atoms
from narupa.imd import ParticleInteraction

from narupatools.ase.constraints.imd_constraint import InteractionConstraint
from narupatools.imd.feature import DynamicsSupportsInteractions, InteractionFeature


class ASEDynamicsSupportsInteractions(DynamicsSupportsInteractions, Protocol):
    """Wrapper around an ASE atoms object."""

    @property
    def atoms(self) -> Atoms:
        """Wrapped ASE Atoms."""
        ...


class ASEIMDFeature(
    InteractionFeature[ASEDynamicsSupportsInteractions, InteractionConstraint]
):
    """Interactive Molecular Dynamics manager for ASE dynamics."""

    @property
    def _system_size(self) -> int:
        return len(self.dynamics.atoms)

    def _make_interaction(
        self, *, key: str, interaction: ParticleInteraction, start_time: float
    ) -> InteractionConstraint:
        constraint = InteractionConstraint(
            dynamics=self._dynamics,
            key=key,
            interaction=interaction,
            start_time=start_time,
        )
        self.dynamics.atoms.constraints.append(constraint)
        return constraint

    def remove_interaction(self, key: str) -> InteractionConstraint:  # noqa: D102
        constraint = super().remove_interaction(key)
        self.dynamics.atoms.constraints.remove(constraint)
        return constraint
