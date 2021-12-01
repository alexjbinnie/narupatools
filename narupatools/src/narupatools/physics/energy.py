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

"""Methods relating to energies and work."""
from typing import Union

import numpy as np

from narupatools.physics.integrate import (
    cumulative_vector_line_integral,
    vector_line_integral,
    vector_line_integral_per_step,
)


def kinetic_energy(
    *, masses: np.ndarray, velocities: np.ndarray
) -> Union[float, np.ndarray]:
    r"""
    Calculate the translational kinetic energy of a system of particles.

    The kinetic energy :math:`K` is defined as:

    .. math:: K = \frac{1}{2} \sum_i m_i | \vec v_i |^2

    :param masses: Masses :math:`m_i` of each particle.
    :param velocities: Velocities :math:`\vec v_i` of each particle.
    """
    return 0.5 * (masses * (velocities ** 2).sum(axis=-1)).sum(axis=-1)  # type: ignore


def total_work(
    *, forces: np.ndarray, positions: np.ndarray, time_axis: int = -2
) -> Union[float, np.ndarray]:
    r"""
    Calculate the total work performed by a force.

    The work is defined by the line integral:

    .. math:: W = \int \vec F \cdot d \vec s

    This is discretized into a sum using the trapezoid rule:

    .. math:: W = \frac{1}{2} \sum_k (\vec F(t_k) + \vec F(t_{k+1}) ) \cdot (\vec r(t_{k+1}) - \vec r(t_{k}) )

    :param forces: Forces :math:`\vec F_i` of each particle.
    :param positions: Velocities :math:`\vec r_i` of each particle.
    :param time_axis: Axis of the forces and positions arrays that represent time.
    """
    return vector_line_integral(forces, positions, 1 + time_axis)


def work_per_step(
    *, forces: np.ndarray, positions: np.ndarray, time_axis: int = -2
) -> np.ndarray:
    r"""
    Calculate the work performed by a force at each timestep.

    The work is defined by the line integral:

    .. math:: W = \int \vec F \cdot d \vec s

    This is discretized into a sum using the trapezoid rule. The work for each step is therefore given by:

    .. math:: W_{k} = \frac{1}{2} (\vec F(t_k) + \vec F(t_{k-1}) ) \cdot (\vec r(t_{k}) - \vec r(t_{k-1}) )

    The first value :math:`W_0` is :math:`0`.

    The axis parameter determines the axis of the forces and position array that represents individual timesteps.
    If the forces and position arrays are of shape [..., timestep, component], then the default axis of -2 should be
    used. If the forces and positions cover a set of particles, and hence are of the form
    [..., timestep, particle, component], then an axis of -3 should be used.

    :param forces: Forces :math:`\vec F_i` of each particle.
    :param positions: Velocities :math:`\vec r_i` of each particle.
    :param time_axis: Axis of the forces and positions arrays that represent time.
    """
    if time_axis > -1:
        raise ValueError("Axis must be less than or equal to -2.")
    return np.insert(  # type: ignore
        vector_line_integral_per_step(forces, positions, axis=1 + time_axis),
        0,
        0,
        axis=1 + time_axis,
    )


def cumulative_work(
    *, forces: np.ndarray, positions: np.ndarray, time_axis: int = -2
) -> np.ndarray:
    r"""
    Calculate the cumulative work performed by a force at each timestep.

    The work is defined by the line integral:

    .. math:: W = \int \vec F \cdot d \vec s

    This is discretized into a sum using the trapezoid rule:

    .. math:: W_{\text{cumulative}}(t_n) = \frac{1}{2} \sum_{k=0}^n (\vec F(t_k) + \vec F(t_{k+1}) ) \cdot (\vec r(t_{k+1}) - \vec r(t_{k}) )

    :param forces: Forces :math:`\vec F_i` of each particle.
    :param positions: Velocities :math:`\vec r_i` of each particle.
    :param time_axis: Axis of the forces and positions arrays that represent time.
    """
    return cumulative_vector_line_integral(forces, positions, axis=1 + time_axis)
