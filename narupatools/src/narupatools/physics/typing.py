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

"""Types for annotating physics methods in a more meaningful way."""

from typing import Sequence, Union

import numpy as np
import numpy.typing as npt

Matrix3x3Like = Union[
    npt.NDArray[np.float64],
    Sequence[npt.NDArray[np.float64]],
    Sequence[Sequence[float]],
]
"""Any type which can be interpreted as a (3, 3) float array by NumPy."""

QuaternionLike = Union[npt.NDArray[np.float64], Sequence[float]]
"""Any type which can be interpreted as a (4,) float array by NumPy."""

Vector3Like = Union[npt.NDArray[np.float64], Sequence[float]]
"""Any type which can be interpreted as a (3,) float array by NumPy."""

VectorNLike = Union[npt.NDArray[np.float64], Sequence[float]]
"""Any type which can be interpreted as a (N,) float array by NumPy."""

ScalarArrayLike = Union[npt.NDArray[np.float64], Sequence[float]]
"""Any type which can be interpreted as a (n,) float array by NumPy."""

Vector3ArrayLike = Union[
    npt.NDArray[np.float64],
    Sequence[npt.NDArray[np.float64]],
    Sequence[Sequence[float]],
]
"""Any type which can be interpreted as a (n, 3) float array by NumPy."""

RotationLike = Union[npt.NDArray[np.float64], Sequence[float]]
"""
Any type which can be interpreted as a rotation.

This could be either a (3, 3) float array or a (4, ) float array.
"""

ScalarArray = npt.NDArray[np.float64]
"""A NumPy float array of shape (n,), representing a list of scalar values."""

IntArray = npt.NDArray[np.int64]
"""A NumPy float array of shape (n,), representing a list of integer values."""

Vector3 = npt.NDArray[np.float64]
"""A NumPy float array of shape (3, ), representing a three dimensional vector."""

Vector3Array = npt.NDArray[np.float64]
"""A NumPy float array of shape (N, 3), representing a list of 3D vectors."""

VectorN = npt.NDArray[np.float64]
"""A NumPy float array of shape (N,), representing an N-dimensional vector."""

Matrix3x3 = npt.NDArray[np.float64]
"""A NumPy float array of shape (3, 3), representing a 3x3 matrix."""

Quaternion = npt.NDArray[np.float64]
"""
A NumPy float array of shape (4,), representing a quaternion.

The quaternion :math:`q = x i + y j + z k + w` is assumed to be in the form as
(x, y, z, w).
"""
