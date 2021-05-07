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

"""An ASE constraint that applies an interactive force."""

from __future__ import annotations

from typing import Optional, Protocol, Tuple, Union

import numpy as np
from ase.atoms import Atoms
from narupa.imd import ParticleInteraction

from narupatools.ase.units import UnitsASE
from narupatools.core import UnitsNarupa
from narupatools.imd import Interaction, calculate_imd_force
from narupatools.physics.typing import Vector3Array
from .constraint import ASEEnergyConstraint
from .constraint import ASEMomentaConstraint

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


class InteractionConstraint(
    Interaction[ASEAtomsWrapper], ASEEnergyConstraint, ASEMomentaConstraint
):
    """
    An ASE constraint that applies an iMD force.

    This conflates both an implementation of a narupatools interaction and an ASE
    constraint.
    """

    def __init__(
        self,
        *,
        dynamics: Union[Atoms, ASEAtomsWrapper],
        key: str,
        interaction: ParticleInteraction,
        start_time: float,
    ):
        """
        Create an ASE constraint that will apply an interactive force.

        :param dynamics: Dynamics this interaction is attached to.
        :param key: Unique key to identify this interaction.
        :param interaction: Initial parameters of the interaction.
        :param start_time: Start time of interaction in simulation time in picoseconds
        """
        if isinstance(dynamics, Atoms):
            dynamics = _ASEAtomsWrapper(dynamics)
        super().__init__(
            dynamics=dynamics, key=key, interaction=interaction, start_time=start_time
        )
        self._atoms: Optional[Atoms] = dynamics.atoms

    def _invalidate_forces_and_energy(self) -> None:
        self._atoms = None
        self._energy = None
        self._forces = None

    def update_energy_and_forces(
        self, *, atoms: Optional[Atoms] = None
    ) -> Tuple[Vector3Array, float]:  # noqa: D102
        if atoms is None:
            atoms = self._dynamics.atoms
        if (
            self._energy is not None
            and self._forces is not None
            and atoms is self._atoms
        ):
            return self._forces, self._energy
        self._atoms = atoms
        return calculate_imd_force(
            interaction=self.interaction,
            positions=atoms.positions[self.particle_indices] * _ASEToNarupa.length,
            masses=atoms.get_masses()[self.particle_indices] * _ASEToNarupa.mass,
        )

    def get_positions(self) -> Vector3Array:  # noqa: D102
        return (  # type:ignore
            self._dynamics.atoms.get_positions()[self.particle_indices]
            * _ASEToNarupa.length
        )

    def adjust_positions(
        self, atoms: Atoms, positions: np.ndarray, /
    ) -> None:  # noqa: D102
        # Assume all interactions depend on positions
        if atoms is self._atoms:
            self._invalidate_forces_and_energy()

    def adjust_momenta(
        self, atoms: Atoms, momenta: np.ndarray, /
    ) -> None:  # noqa: D102
        # Assume no interactions depend on velocities
        # When they do, this should conditionally invalidate the cache
        pass

    def adjust_forces(self, atoms: Atoms, forces: np.ndarray, /) -> None:  # noqa: D102
        self._forces, self._energy = self.update_energy_and_forces(atoms=atoms)
        forces[self._interaction.particles] += self.forces * _NarupaToASE.force

    def adjust_potential_energy(self, /, atoms: Atoms) -> float:  # noqa: D102
        self._forces, self._energy = self.update_energy_and_forces(atoms=atoms)
        return self.potential_energy * _NarupaToASE.energy
