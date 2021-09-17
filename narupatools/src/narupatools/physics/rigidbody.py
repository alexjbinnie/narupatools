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

"""
Methods for dealing with systems of particles, such as center of mass.

Methods here are not specific with units - as long as arguments are provided in
consistent units, then the calculated result will be correct.
"""
import math
from typing import Optional

import numpy as np
from numpy.linalg import inv

from .matrix import zero_matrix
from .typing import (
    Matrix3x3,
    ScalarArray,
    ScalarArrayLike,
    Vector3,
    Vector3Array,
    Vector3ArrayLike,
    Vector3Like,
)
from .vector import (
    cross_product,
    left_vector_triple_product_matrix,
    sqr_magnitude,
    zero_vector,
)


def center_of_mass(*, masses: ScalarArrayLike, positions: Vector3ArrayLike) -> Vector3:
    r"""
    Calculate the center of mass :math:`R` of a collection of particles.

    This is defined as the mass weighted average position:

    .. math:: R = \frac{\sum_i m_i r_i}{\sum_i m_i}

    :param masses: List of masses :math:`m_i` of each particle.
    :param positions: List of positions :math:`r_i` of each particle.
    :raises ValueError: Masses and position arrays are different lengths.
    :return: Center of mass of the particles, in the same units as positions was
             provided in.
    """
    positions: Vector3Array = np.asfarray(positions)
    masses: ScalarArray = np.asfarray(masses)

    count = len(positions)
    if len(masses) != count:
        raise ValueError(
            f"Mismatch between number of positions ({len(positions)}) "
            f"and number of masses ({len(masses)})"
        )

    total_center = zero_vector()
    for i in range(count):
        total_center += masses[i] * positions[i]
    return total_center / sum(masses)


def center_of_mass_velocity(
    *, masses: ScalarArrayLike, velocities: Vector3ArrayLike
) -> Vector3:
    r"""
    Calculate the velocity of the center of mass of a collection of particles.

    This is defined as the mass weighted average velocity:

    .. math:: V_R = \frac{\sum_i m_i v_i}{\sum_i v_i}

    :param masses: List of masses :math:`m_i` of each particle.
    :param velocities: List of velocities :math:`v_i` of each particle.
    :return: Velocity of the center of mass of the particles, in units of [distance]
             / [time].
    """
    # Reuse center of mass calculation, as its the same with positions exchanged for
    # velocities
    return center_of_mass(masses=masses, positions=velocities)


def center_of_mass_acceleration(
    *, masses: ScalarArrayLike, accelerations: Vector3ArrayLike
) -> Vector3:
    r"""
    Calculate the acceleration of the center of mass of a collection of particles.

    This is defined as the mass weighted average acceleration:

    .. math:: A_R = \frac{\sum_i m_i a_i}{\sum_i a_i}

    :param masses: List of masses :math:`m_i` of each particle.
    :param accelerations: List of accelerations :math:`a_i` of each particle.
    :return: Acceleration of the center of mass of the particles, in units of [distance]
             / [time] ** 2.
    """
    # Reuse center of mass calculation, as its the same with positions exchanged for
    # accelerations
    return center_of_mass(masses=masses, positions=accelerations)


def spin_angular_momentum(
    *, masses: ScalarArray, positions: Vector3Array, velocities: Vector3Array
) -> Vector3:
    r"""
    Calculate the spin angular momentum of a set of particles.

    This is defined as:

    .. math:: L = \sum_i m_i r_i \times v_i

    where :math:`r_i` and :math:`v_i` are the particles' position and velocity relative
    to the center of mass.

    :param masses: List of masses :math:`m_i` of each particle.
    :param positions: List of positions :math:`R_i` of each particle.
    :param velocities: List of velocities :math:`V_i` of each particle.
    :return: Spin angular momentum :math:`L` in units of [mass] * [distance] squared /
             [time].
    """
    com = center_of_mass(masses=masses, positions=positions)
    com_velocity = center_of_mass_velocity(masses=masses, velocities=velocities)
    angular_momentum = np.array([0.0, 0.0, 0.0], dtype=float)
    for i in range(len(masses)):
        angular_momentum += masses[i] * cross_product(
            positions[i] - com, velocities[i] - com_velocity
        )
    return angular_momentum


def orbital_angular_momentum(
    *,
    masses: ScalarArray,
    positions: Vector3Array,
    velocities: Vector3Array,
    origin: Optional[Vector3Like] = None,
) -> Vector3:
    r"""
    Calculate the orbital angular momentum of a set of particles about an origin.

    This is defined as:

    .. math:: L = M (R - c) \times V

    where :math:`M`, :math:`R` and :math:`V` are the total mass, center of mass and the
    velocity of the center of mass respectively.

    The origin defaults to the origin (0, 0, 0).

    :param masses: List of masses :math:`m_i` of each particle.
    :param positions: List of positions :math:`R_i` of each particle.
    :param velocities: List of velocities :math:`V_i` of each particle.
    :param origin: Origin about where to calculate the orbital angular momentum.
    :return: Orbital angular momentum :math:`L` in units of [mass] * [distance] squared
             / [time].
    """
    if origin is None:
        origin = zero_vector()
    else:
        origin = np.asfarray(origin)
    com = center_of_mass(masses=masses, positions=positions)
    com_velocity = center_of_mass_velocity(masses=masses, velocities=velocities)
    total_mass = np.sum(np.asfarray(masses))
    return total_mass * np.cross(com - origin, com_velocity)  # type: ignore


def radius_of_gyration(
    *,
    masses: ScalarArray,
    positions: Vector3Array,
    axis: Vector3,
    origin: Optional[Vector3Like] = None,
) -> float:
    r"""
    Calculate the radius of gyration about an axis.

    The radius of gyration :math:`R` about an axis is the radius at which a fictitious particle of total mass :math:`M`
    would have the same moment of inertia as the provided system does about the axis and origin.

    It is the solution to the equation:

    .. math:: M R^2 = I

    where :math:`M` is the total mass and :math:`I` is the moment of inertia about the axis (see
    :func:`moment_of_inertia`) for more details.

    By default, the origin is taken to be the center of mass.

    :param masses: List of masses :math:`m_i` of each particle.
    :param positions: List of positions :math:`R_i` of each particle.
    :param axis: Axis about which to calculate the radius of gyration.
    :param origin: Optional origin to calculate the moment of inertia around. Defaults
                   to the center of mass.
    :return: Radius of gyration :math:`R` with respect to the given axis and origin, in
             units of [distance].
    """
    return math.sqrt(
        moment_of_inertia(masses=masses, positions=positions, axis=axis, origin=origin)
        / masses.sum()
    )


def moment_of_inertia(
    *,
    masses: ScalarArray,
    positions: Vector3Array,
    axis: Vector3,
    origin: Optional[Vector3Like] = None,
) -> float:
    r"""
    Calculate the moment of inertia of a set of particles about an axis.

    The moment of inertia acts similarly to a rotational analogue to mass, and describes the distribution of mass
    about the provided axis.

    .. math:: I = \sum_i m_i |r_i^\perp|^2

    where :math:`m_i` is the mass of each particle and :math:`r_i^\perp` is the distance of the particle from the axis.

    By default, the origin is taken to be the center of mass.

    :param masses: List of masses :math:`m_i` of each particle.
    :param positions: List of positions :math:`R_i` of each particle.
    :param axis: Axis about which to calculate the radius of gyration.
    :param origin: Optional origin to calculate the moment of inertia around. Defaults
                   to the center of mass.
    :return: Moment of inertia :math:`I` with respect to the given axis and origin, in
             units of [mass] * [distance] squared
    """
    tensor = moment_of_inertia_tensor(masses=masses, positions=positions, origin=origin)
    return np.dot(axis, tensor @ axis) / np.dot(axis, axis)  # type: ignore[no-any-return]


def moment_of_inertia_tensor(
    *,
    masses: ScalarArray,
    positions: Vector3Array,
    origin: Optional[Vector3Like] = None,
) -> Matrix3x3:
    r"""
    Calculate the moment of inertia tensor of a set of particles.

    This matrix fulfills the identity:

    .. math:: I v = \sum_i m_i r_i \times (v \times r_i)

    for all vectors :math:`v`. Here, :math:`r_i` is the position of the i-th particle
    relative to the center of mass.

    By default, the origin is taken to be the center of mass.

    :param masses: List of masses :math:`m_i` of each particle.
    :param positions: List of positions :math:`R_i` of each particle.
    :param origin: Optional origin to calculate the moment of inertia around. Defaults
                   to the center of mass.
    :return: Moment of inertia tensor :math:`I` with respect to the center of mass, in
             units of [mass] * [distance] squared
    """
    if origin is None:
        origin = center_of_mass(masses=masses, positions=positions)
    else:
        origin = np.asfarray(origin)
    tensor = zero_matrix()
    for i in range(len(masses)):
        tensor -= masses[i] * left_vector_triple_product_matrix(
            positions[i] - origin, positions[i] - origin
        )
    return tensor


def angular_velocity(
    *,
    masses: ScalarArray,
    positions: Vector3Array,
    velocities: Vector3Array,
) -> Vector3:
    r"""
    Calculate the angular velocity of a set of particles.

    This is performed by calculating the angular momentum about the center of mass, and
    inverting the inertia tensor to obtain :math:`\omega`:

    .. math:: \omega = I^{-1} L

    :param masses: List of masses :math:`m_i` of each particle.
    :param positions: List of positions :math:`R_i` of each particle.
    :param velocities: List of velocities :math:`V_i` of each particle.
    :return: Angular velocity :math:`\omega` in units of radians / [time].
    """
    L = spin_angular_momentum(masses=masses, positions=positions, velocities=velocities)
    inertia = moment_of_inertia_tensor(masses=masses, positions=positions)
    return inv(inertia) @ L  # type: ignore


def distribute_angular_velocity(
    *,
    masses: Optional[ScalarArrayLike] = None,
    positions: Vector3ArrayLike,
    angular_velocity: Vector3Like,
    origin: Optional[Vector3Like] = None,
) -> Vector3Array:
    r"""
    Calculate the velocities that correspond to an angular velocity about an origin.

    This is defined by:

    .. math:: v_i = \omega \times (r_i - c)

    The origin :math:`c` defaults to the center of mass of the particles.

    :param masses: List of masses :math:`m_i` of each particle.
    :param positions: List of positions :math:`R_i` of each particle.
    :param angular_velocity: Angular velocity :math:`\omega` to assign to the system.
    :param origin: Origin :math:`c` to calculate velocity about.
    :raises ValueError: Neither an origin or masses are provided.
    :return: Velocities :math:`V_i` in units of [distance] / [time].
    """
    if origin is None:
        if masses is None:
            raise ValueError("Either an origin or masses must be provided.")
        o = center_of_mass(masses=masses, positions=positions)
    else:
        o = np.asfarray(origin)
    velocities = [np.cross(angular_velocity, position - o) for position in positions]
    return np.asfarray(velocities)  # type: ignore


def kinetic_energy(*, masses: ScalarArrayLike, velocities: Vector3ArrayLike) -> float:
    r"""
    Calculate the total kinetic energy of a set of particles.

    This is given by:

    .. math:: K = \frac{1}{2} \sum_i m_i v_i^2

    :param masses: List of masses :math:`m_i` of each particle.
    :param velocities: List of velocities :math:`v_i` of each particle.
    :return: Kinetic energy :math:`K` in units of [distance] / [time].
    """
    return sum(
        0.5 * mass * sqr_magnitude(velocity)
        for mass, velocity in zip(masses, velocities)
    )
