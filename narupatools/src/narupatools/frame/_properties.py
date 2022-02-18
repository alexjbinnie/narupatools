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

from typing import Any, Collection, Iterator, Protocol, Union

import numpy as np
import numpy.typing as npt
from scipy.spatial.transform import Rotation

from narupatools.override import override
from narupatools.physics import quaternion
from narupatools.physics.rigidbody import (
    angular_velocity,
    center_of_mass,
    center_of_mass_velocity,
    moment_of_inertia_tensor,
    principal_axes,
)
from narupatools.physics.thermodynamics import maxwell_boltzmann_velocities
from narupatools.physics.typing import (
    ScalarArray,
    Vector3,
    Vector3Array,
    Vector3ArrayLike,
)

from ._select import select


class StaticStructureProperties(Protocol):
    """Base protocol for an object which represents a static system of particles."""

    @property
    def positions(self) -> Vector3Array:
        """Positions of particles in nanometers."""
        raise AttributeError

    @property
    def masses(self) -> ScalarArray:
        """Masses of particles in daltons."""
        raise AttributeError

    @property
    def orientations(self) -> npt.NDArray[quaternion]:
        """Orientations of each atom as unit quaternions."""
        raise AttributeError

    @property
    def moments_of_inertia(self) -> Vector3Array:
        """
        Moments of inertia for each particle abouts its origin in its local frame.

        :return: Array of moments of inertia, either scalars (for symmetric shapes) or
                 3-vectors.
        """
        raise AttributeError


class WritableStaticStructureProperties(StaticStructureProperties, Protocol):
    """Base protocol for an object which represents a static system of particles."""

    @override(StaticStructureProperties.positions)
    @property
    def positions(self) -> Vector3Array:  # noqa: D102
        raise AttributeError

    @positions.setter
    def positions(self, value: Vector3ArrayLike) -> None:
        ...


class DynamicStructureMethods:
    """Mixin methods for calculating properties."""

    def center_of_mass_velocity(self: DynamicStructureProperties) -> Vector3:
        """Calculate the velocity of the center of mass."""
        return center_of_mass_velocity(velocities=self.velocities, masses=self.masses)

    def angular_velocity(self: DynamicStructureProperties) -> Vector3:
        """Calculate the angular velocity about the center of mass."""
        return angular_velocity(
            positions=self.positions, velocities=self.velocities, masses=self.masses
        )

    def select(
        self: DynamicStructureProperties, selection: Union[str, np.ndarray]
    ) -> SelectionView:
        """Select a subset using either indices or selections."""
        return SelectionView(self, selection)

    def center_of_mass(self: StaticStructureProperties) -> Vector3:
        """Calculate the center of mass."""
        return center_of_mass(positions=self.positions, masses=self.masses)

    def moment_of_inertia_tensor(self: StaticStructureProperties) -> np.ndarray:
        """Calculate the moment of inertia tensor."""
        return moment_of_inertia_tensor(positions=self.positions, masses=self.masses)

    def principal_axes(self: StaticStructureProperties) -> Vector3Array:
        """Calculate the principal axes."""
        return principal_axes(positions=self.positions, masses=self.masses)

    def translate_to(
        self: WritableStaticStructureProperties, position: Vector3
    ) -> None:
        """Translate the system so the center of mass is shifted to the provided position."""
        com = self.center_of_mass()  # type: ignore
        self.positions += position - com

    def randomly_orient(self: WritableStaticStructureProperties) -> None:
        """Randomly orientate the object around its center of mass."""
        com = self.center_of_mass()  # type: ignore
        rot = Rotation.random()
        self.positions = com + rot.apply(self.positions - com)

    def distribute_velocities(
        self: WritableDynamicStructureProperties, temperature: float
    ) -> None:
        """
        Distribute velocities using the Maxwell-Boltzmann distribution.

        :param temperature: Temperature of the Maxwell-Boltzmann distribution in Kelvin.
        """
        self.velocities = maxwell_boltzmann_velocities(
            masses=self.masses, temperature=temperature
        )


class DynamicStructureProperties(StaticStructureProperties, Protocol):
    """Base protocol for an object which represents a dynamic system of particles."""

    @property
    def velocities(self) -> Vector3Array:
        """Velocities of particles in nanometers per picosecond."""
        raise AttributeError

    @property
    def momenta(self) -> Vector3Array:
        """Momenta of particles in dalton nanometers per picosecond."""
        return self.velocities * self.masses

    @property
    def forces(self) -> Vector3Array:
        """Forces on particles in kilojoules per mole per nanometer."""
        raise AttributeError

    @property
    def kinetic_energy(self) -> float:
        """Kinetic energy in kilojoules per mole."""
        raise AttributeError

    @property
    def potential_energy(self) -> float:
        """Potential energy in kilojoules per mole."""
        raise AttributeError

    @property
    def angular_momenta(self) -> Vector3Array:
        """
        Angular momentum of each particle abouts its center of mass.

        :raises AttributeError: Angular momenta is not defined for this system.
        :return: Array of angular momenta in dalton nanometer squared per picosecond.
        """
        raise AttributeError

    @property
    def angular_velocities(self) -> Vector3Array:
        """
        Angular velocity of each particle abouts its center of mass.

        :return: Array of angular velocities in radians per picoseconds.
        """
        raise AttributeError

    @property
    def torques(self) -> Vector3Array:
        """
        Torques on each particle abouts its center of mass.

        :return: Array of torques in kilojoules per mole.
        """
        raise AttributeError


class WritableDynamicStructureProperties(DynamicStructureProperties, Protocol):
    """Base protocol for an object which represents a dynamic system of particles."""

    @override(DynamicStructureProperties.velocities)
    @property
    def velocities(self) -> Vector3Array:  # noqa: D102
        raise AttributeError

    @velocities.setter
    def velocities(self, value: Vector3ArrayLike) -> None:
        ...


class SelectionView(
    DynamicStructureProperties, DynamicStructureMethods, Collection[float]
):
    """
    View to a subset of a system, such as a FrameData or dynamics.

    A selection consists of a set of indices, representing the indices of the particles within the
    selection, and an object that this selection is chosen from.

    Selections from the same object may be combined using standard set operations, such as &, | and
    ^. They may also be compared to each other using set notation for subsets and supersets.
    """

    def __init__(
        self, source: DynamicStructureProperties, selection: Union[str, np.ndarray]
    ):
        self._source = source
        if isinstance(selection, str):
            selection = select(source, selection)
        self._selection = selection

    def __iter__(self) -> Iterator[float]:
        return iter(self._selection)

    def __contains__(self, item: Any) -> bool:
        return item in self._selection

    def __len__(self) -> int:
        return len(self._selection)

    def __le__(self, other: Any) -> bool:
        return set(self._selection) <= set(other)

    def __lt__(self, other: Any) -> bool:
        return set(self._selection) < set(other)

    def __eq__(self, other: Any) -> bool:
        return set(self._selection) == set(other)

    def __ne__(self, other: Any) -> bool:
        return set(self._selection) != set(other)

    def __gt__(self, other: Any) -> bool:
        return set(self._selection) > set(other)

    def __ge__(self, other: Any) -> bool:
        return set(self._selection) >= set(other)

    def __and__(self, other: Any) -> SelectionView:
        if not isinstance(other, SelectionView):
            return NotImplemented
        if self._source != other._source:
            raise ValueError("Cannot intersect two selections from different objects.")
        return SelectionView(
            self._source, np.intersect1d(self._selection, other._selection)
        )

    def __or__(self, other: Any) -> SelectionView:
        if not isinstance(other, SelectionView):
            return NotImplemented
        if self._source != other._source:
            raise ValueError(
                "Cannot take the union of two selections from different objects."
            )
        return SelectionView(
            self._source, np.union1d(self._selection, other._selection)
        )

    def __xor__(self, other: Any) -> SelectionView:
        if not isinstance(other, SelectionView):
            return NotImplemented
        if self._source != other._source:
            raise ValueError(
                "Cannot take the symmetric difference of two selections from different objects."
            )
        return SelectionView(
            self._source, np.setxor1d(self._selection, other._selection)
        )

    def __sub__(self, other: Any) -> SelectionView:
        if not isinstance(other, SelectionView):
            return NotImplemented
        if self._source != other._source:
            raise ValueError(
                "Cannot take the set difference of two selections from different objects."
            )
        return SelectionView(
            self._source, np.setdiff1d(self._selection, other._selection)
        )

    @property
    def indices(self) -> np.ndarray:
        """Indices of the particles that belong in this selection."""
        return self._selection.copy()

    @override(DynamicStructureProperties.masses)
    @property
    def masses(self) -> ScalarArray:  # noqa: D102
        return self._source.masses[..., self._selection]

    @masses.setter
    def masses(self, value: ScalarArray) -> None:
        masses = self._source.masses
        masses[..., self._selection] = value
        self._source.masses = masses  # type: ignore

    @override(DynamicStructureProperties.positions)
    @property
    def positions(self) -> Vector3Array:  # noqa: D102
        return self._source.positions[..., self._selection, :]

    @positions.setter
    def positions(self, value: Vector3Array) -> None:
        positions = self._source.positions
        positions[..., self._selection, :] = value
        self._source.positions = positions  # type: ignore

    @override(DynamicStructureProperties.velocities)  # type: ignore
    @property
    def velocities(self) -> Vector3Array:  # type: ignore  # noqa: D102
        return self._source.velocities[..., self._selection, :]

    @velocities.setter
    def velocities(self, value: Vector3Array) -> None:
        velocities = self._source.velocities
        velocities[..., self._selection, :] = value
        self._source.velocities = velocities  # type: ignore

    @override(DynamicStructureProperties.orientations)
    @property
    def orientations(self) -> npt.NDArray[quaternion]:  # noqa: D102
        return self._source.orientations[..., self._selection, :]

    @override(DynamicStructureProperties.moments_of_inertia)
    @property
    def moments_of_inertia(self) -> Vector3Array:  # noqa: D102
        """
        Moments of inertia for each particle abouts its origin in its local frame.

        :return: Array of moments of inertia, either scalars (for symmetric shapes) or
                 3-vectors.
        """
        raise AttributeError

    def __repr__(self) -> str:
        return f'<SelectionView selection="{self._selection}" of {self._source}>'

    def __array__(self, dtype: npt.DTypeLike = None) -> np.ndarray:
        if dtype is None:
            return self.indices.copy()
        else:
            return np.asarray(self.indices, dtype=dtype)
