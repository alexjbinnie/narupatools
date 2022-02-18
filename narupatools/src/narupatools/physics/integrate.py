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

"""Calculate numeric integrals and line integrals."""

from typing import Any, List, Union

import numpy as np


def integral(
    integrand: np.ndarray, variable: np.ndarray, /, axis: int = -1
) -> Union[float, np.ndarray]:
    r"""
    Integrate a function :math:`f(x)` over a variable :math:`x`.

    This function acts on a set of samples, with arrays of :math:`f(x)` and :math:`x`. The
    integral is evaluated by the trapezoid rule:

    .. math:: \int f(x) \ dx \approx \frac12 \sum_i (f(x_i) + f(x_{i+1})) (x_{i+1} - x_i)

    :param integrand: Values of the function to be integrated.
    :param variable: Corresponding integration variable at the same samples.
    :return: Value of the integration.
    """
    integrand = np.asfarray(integrand)
    variable = np.asfarray(variable)

    # How many axes exist between the samples column and the final column (the components)
    # The default (-1) indicates there is only one axis (the components)
    end_indices = -1 - axis
    last_axes: List[slice] = [np.s_[:]] * end_indices
    # Index selection for the starts and ends of each segment
    step_start: Any = (..., np.s_[:-1], *last_axes)
    step_end: Any = (..., np.s_[1:], *last_axes)

    return 0.5 * (  # type: ignore
        (integrand[step_start] + integrand[step_end])
        * (variable[step_end] - variable[step_start])
    ).sum(axis=axis)


def vector_line_integral(
    integrand: np.ndarray, variable: np.ndarray, /, axis: int = -1
) -> Union[float, np.ndarray]:
    r"""
    Integrate a function :math:`\vec f(\vec x)` along a line integral in :math:`\vec x`.

    This function acts on a set of samples, with arrays of :math:`\vec f(\vec x)` and :math:`\vec x`. The
    integral is evaluated by the trapezoid rule:

    .. math:: \int \vec f(\vec x) \cdot d \vec x \approx \frac12 \sum_i (\vec f(\vec x_i) + \vec f(\vec x_{i+1})) \cdot (\vec x_{i+1} - \vec x_i)

    The default behaviour with the parameter axis equal to -1 is to consider the integrand and variable arrays to have shape:

    * integrand: (..., samples, components)
    * variable: (..., samples, components)

    The axis parameter of -1 corresponds to the samples column, as the components column has already been elimited by the dot product.

    The axis parameter should be set to -2 if there is some level of grouping column after the samples column, for example:

    * integrand: (..., samples, group, components)
    * variable: (..., samples, group, components)

    :param integrand: Values of the function to be integrated.
    :param variable: Corresponding integration variable at the same samples.
    :param axis: Axis of the integrand and variable arrays that represents the discrete samples, ignoring the last
                 array.
    :return: Value of the integration.
    """
    return vector_line_integral_per_step(integrand, variable, axis=axis).sum(axis=axis)  # type: ignore


def cumulative_vector_line_integral(
    integrand: np.ndarray, variable: np.ndarray, /, axis: int = -1
) -> np.ndarray:
    r"""
    Integrate a function :math:`\vec f(\vec x)` along a line integral in :math:`\vec x`, calculating the cumulative sum.

    This function acts on a set of samples, with arrays of :math:`\vec f(\vec x)` and :math:`\vec x`. The
    integral is evaluated by the trapezoid rule:

    .. math:: \int \vec f(\vec x) \cdot d \vec x \approx \frac12 \sum_i (\vec f(\vec x_i) + \vec f(\vec x_{i+1})) \cdot (\vec x_{i+1} - \vec x_i)

    :param integrand: Values of the function to be integrated.
    :param variable: Corresponding integration variable at the same samples.
    :param axis: Axis of the integrand and variable arrays that represents the discrete samples, ignoring the last
                 array.
    :return: Value of the integration.
    """
    cumm = np.cumsum(
        vector_line_integral_per_step(integrand, variable, axis=axis), axis=axis
    )
    return np.insert(cumm, 0, 0, axis)


def vector_line_integral_per_step(
    integrand: np.ndarray, variable: np.ndarray, /, axis: int = -1
) -> np.ndarray:
    r"""
    Calculate the contribution of each step to a line integral of :math:`\vec f(\vec x(t))` along a curve :math:`\vec x(t)`.

    This function acts on a set of samples, with arrays of :math:`\vec f(\vec x)` and :math:`\vec x`. The
    integral is evaluated by the trapezoid rule:

    .. math:: \int \vec f(\vec x) \cdot d \vec x \approx \frac12 \sum_i (\vec f(\vec x_i) + \vec f(\vec x_{i+1})) \cdot (\vec x_{i+1} - \vec x_i)

    The default behaviour with the parameter axis equal to -1 is to consider the integrand and variable arrays to have shape:

    * integrand: (..., samples, components)
    * variable: (..., samples, components)

    The axis parameter of -1 corresponds to the samples column, as the components column has already been elimited by the dot product.

    The axis parameter should be set to -2 if there is some level of grouping column after the samples column, for example:

    * integrand: (..., samples, group, components)
    * variable: (..., samples, group, components)

    :param integrand: Values of the function to be integrated.
    :param variable: Corresponding integration variable at the same samples.
    :param axis: Axis of the integrand and variable arrays that represents the discrete samples, ignoring the last
                 array.
    :return: Value of the integration.
    """
    integrand = np.asfarray(integrand)
    variable = np.asfarray(variable)

    # How many axes exist between the samples column and the final column (the components)
    # The default (-1) indicates there is only one axis (the components)
    end_indices = -axis
    last_axes: List[slice] = [np.s_[:]] * end_indices
    # Index selection for the starts and ends of each segment
    step_start: Any = (..., np.s_[:-1], *last_axes)
    step_end: Any = (..., np.s_[1:], *last_axes)

    return 0.5 * (  # type: ignore
        (integrand[step_start] + integrand[step_end])
        * (variable[step_end] - variable[step_start])
    ).sum(axis=-1)
