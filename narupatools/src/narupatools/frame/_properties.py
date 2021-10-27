from __future__ import annotations

from typing import Protocol, Union

import numpy as np
import numpy.typing as npt
from scipy.spatial.transform import Rotation

from narupatools.physics import quaternion
from narupatools.physics.rigidbody import (
    angular_velocity,
    center_of_mass,
    center_of_mass_velocity,
    moment_of_inertia_tensor,
    principal_axes,
)
from narupatools.physics.typing import (
    ScalarArray,
    Vector3,
    Vector3Array,
    Vector3ArrayLike,
)

from ..physics.thermodynamics import maxwell_boltzmann_velocities
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

    @property
    def positions(self) -> Vector3Array:
        raise AttributeError

    @positions.setter
    def positions(self, value: Vector3ArrayLike) -> None:
        ...


class DynamicStructureMethods:
    def center_of_mass_velocity(self: DynamicStructureProperties) -> Vector3:
        return center_of_mass_velocity(velocities=self.velocities, masses=self.masses)

    def angular_velocity(self: DynamicStructureProperties) -> Vector3:
        return angular_velocity(
            positions=self.positions, velocities=self.velocities, masses=self.masses
        )

    def select(
        self: DynamicStructureProperties, selection: Union[str, np.ndarray]
    ) -> SelectionView:
        return SelectionView(self, selection)

    def center_of_mass(self: StaticStructureProperties) -> Vector3:
        return center_of_mass(positions=self.positions, masses=self.masses)

    def moment_of_inertia_tensor(self: StaticStructureProperties) -> np.ndarray:
        return moment_of_inertia_tensor(positions=self.positions, masses=self.masses)

    def principal_axes(self: StaticStructureProperties) -> Vector3Array:
        return principal_axes(positions=self.positions, masses=self.masses)

    def translate_to(
        self: WritableStaticStructureProperties, position: Vector3
    ) -> None:
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


class WritableDynamicStructureProperties(StaticStructureProperties, Protocol):
    """Base protocol for an object which represents a dynamic system of particles."""

    @property
    def velocities(self) -> Vector3Array:
        raise AttributeError

    @velocities.setter
    def velocities(self, value: Vector3ArrayLike) -> None:
        ...


class SelectionView(DynamicStructureProperties, DynamicStructureMethods):
    """View to a subset of a system, such as a FrameData or dynamics."""

    def __init__(
        self, source: DynamicStructureProperties, selection: Union[str, np.ndarray]
    ):
        self._source = source
        if isinstance(selection, str):
            selection = select(source, selection)
        self._selection = selection

    @property
    def indices(self) -> np.ndarray:
        return self._selection.copy()

    @property
    def masses(self) -> ScalarArray:
        return self._source.masses[self._selection]  # type: ignore

    @masses.setter
    def masses(self, value: ScalarArray) -> None:
        masses = self._source.masses
        masses[self._selection] = value
        self._source.masses = masses  # type: ignore

    @property
    def positions(self) -> Vector3Array:
        return self._source.positions[self._selection]  # type: ignore

    @positions.setter
    def positions(self, value: Vector3Array) -> None:
        positions = self._source.positions
        positions[self._selection] = value
        self._source.positions = positions  # type: ignore

    @property  # type: ignore
    def velocities(self) -> Vector3Array:  # type: ignore
        return self._source.velocities[self._selection]  # type: ignore

    @velocities.setter
    def velocities(self, value: Vector3Array) -> None:
        velocities = self._source.velocities
        velocities[self._selection] = value
        self._source.velocities = velocities  # type: ignore

    @property
    def orientations(self) -> npt.NDArray[quaternion]:
        return self._source.orientations[self._selection]  # type: ignore

    @property
    def moments_of_inertia(self) -> Vector3Array:
        """
        Moments of inertia for each particle abouts its origin in its local frame.

        :return: Array of moments of inertia, either scalars (for symmetric shapes) or
                 3-vectors.
        """
        raise AttributeError

    def __repr__(self) -> str:
        return f'<SelectionView selection="{self._selection}" of {self._source}>'

    def __array__(self, dtype: npt.DTypeLike = None) -> np.ndarray:
        if dtype == None:
            return self.indices.copy()
        else:
            return np.asarray(self.indices, dtype=dtype)
