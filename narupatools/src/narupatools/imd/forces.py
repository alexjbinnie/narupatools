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

"""Code for calculating IMD forces."""

from typing import Any, Dict, Protocol, Tuple

import numpy as np
from narupa.imd import ParticleInteraction

from narupatools.physics.force import (
    gaussian_force_and_energy,
    mass_weighted_forces,
    spring_force_and_energy,
)
from narupatools.physics.rigidbody import center_of_mass
from narupatools.physics.typing import ScalarArray, Vector3, Vector3Array


def calculate_imd_force(
    *, interaction: ParticleInteraction, positions: Vector3Array, masses: Vector3Array
) -> Tuple[Vector3Array, float]:
    """
    Calculate the forces and energy to apply for a given interaction.

    :param interaction: Interaction to calculate forces and energy for.
    :param positions: Positions of particles affected by this interaction.
    :param masses: Masses of particles affected by this interaction.
    """
    kwargs = {"position": interaction.position, **interaction.properties}

    return _calculate_imd_force(
        positions=positions,
        masses=masses,
        interaction_type=interaction.interaction_type,
        interaction_scale=interaction.scale,
        **kwargs,
    )


class UnknownInteractiveForceError(ValueError):
    """Error raised when an interaction with an unknown type is used."""

    def __init__(self, type: str, /):
        super().__init__(f"Unknown interactive force type: {type}")


class InteractionTypeAlreadyDefinedError(ValueError):
    """Error raised when an interaction type already exists."""

    def __init__(self, type: str, /):
        super().__init__(f"Interaction type {type} already exists.")


class InteractiveForce(Protocol):
    """Protocol for a force function that returns a set of forces and an energy."""

    def __call__(self, **kwargs: Any) -> Tuple[Vector3Array, float]:
        """
        Calculate the per-particle forces and potential energy for this interaction.

        The forces must be in kilojoules per mole per nanometer, and the energy in
        kilojoules per mole.

        :param kwargs: Arbitrary arguments to pass to the force.
        """
        ...


def _constant_force(
    *,
    positions: Vector3Array,
    masses: ScalarArray,
    interaction_scale: float,
    force: Vector3,
    **kwargs: Any,
) -> Tuple[Vector3Array, float]:
    center = center_of_mass(positions=positions, masses=masses)
    force = np.array(force, dtype=float) * interaction_scale
    energy = -np.dot(center, force) * interaction_scale
    return mass_weighted_forces(force=force, masses=masses), energy


def _gaussian_force(
    *,
    positions: Vector3Array,
    masses: ScalarArray,
    interaction_scale: float,
    position: Vector3,
    **kwargs: Any,
) -> Tuple[Vector3Array, float]:
    center = center_of_mass(positions=positions, masses=masses)
    force, energy = gaussian_force_and_energy(offset=center - position)
    force *= interaction_scale
    energy *= interaction_scale
    return mass_weighted_forces(force=force, masses=masses), energy


def _spring_force(
    *,
    positions: Vector3Array,
    masses: ScalarArray,
    interaction_scale: float,
    position: Vector3,
    **kwargs: Any,
) -> Tuple[Vector3Array, float]:
    center = center_of_mass(positions=positions, masses=masses)
    force, energy = spring_force_and_energy(offset=center - position)
    force *= interaction_scale
    energy *= interaction_scale
    return mass_weighted_forces(force=force, masses=masses), energy


CONSTANT_INTERACTION_TYPE = "constant"
SPRING_INTERACTION_TYPE = "spring"
GAUSSIAN_INTERACTION_TYPE = "gaussian"

_GLOBAL_FORCE_FUNCTIONS: Dict[str, InteractiveForce] = {}


def register_imd_force_type(type: str, force: InteractiveForce, /) -> None:
    """
    Register a new kind of interactive force.

    :param type: String key that defines this interaction type.
    :param force: Function that calculates the force and energy.
    :raises InteractionTypeAlreadyDefinedError: Interaction with the given name already
                                                exists.
    """
    if type in _GLOBAL_FORCE_FUNCTIONS:
        raise InteractionTypeAlreadyDefinedError(type)
    _GLOBAL_FORCE_FUNCTIONS[type] = force


register_imd_force_type(GAUSSIAN_INTERACTION_TYPE, _gaussian_force)  # type: ignore
register_imd_force_type(SPRING_INTERACTION_TYPE, _spring_force)  # type: ignore
register_imd_force_type(CONSTANT_INTERACTION_TYPE, _constant_force)  # type: ignore


def get_force_function(type: str, /) -> InteractiveForce:
    """
    Get the force function defined for a given interaction type.

    :param type: Interaction type to look up.
    :raises UnknownInteractiveForceError: Given interaction type is not recognized.
    """
    if type in _GLOBAL_FORCE_FUNCTIONS:
        return _GLOBAL_FORCE_FUNCTIONS[type]
    raise UnknownInteractiveForceError(type)


def _calculate_imd_force(
    *,
    positions: Vector3Array,
    masses: ScalarArray,
    interaction_type: str,
    interaction_scale: float,
    **kwargs: Any,
) -> Tuple[Vector3Array, float]:
    """
    Calculate the forces and energy to apply for a given interaction.

    :param positions: Positions of particles affected by this interaction.
    :param masses: Masses of particles affected by this interaction.
    :param interaction_type: Type of the interaction.
    :param interaction_scale: Scaling factor of the interaction.
    :param kwargs: Arbitrary arguments to pass to the force function.
    :return: Tuple of the forces and energy, in matching units to the units of the
             positions and masses.
    """
    force_function = get_force_function(interaction_type)
    return force_function(
        positions=positions,
        masses=masses,
        interaction_scale=interaction_scale,
        **kwargs,
    )
