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

"""Utility methods for using matrices."""

import numpy as np

from .typing import Matrix3x3, Matrix3x3Like


def zero_matrix() -> Matrix3x3:
    """Create a 3x3 matrix with all zero entries."""
    return np.zeros((3, 3))


def identity_matrix() -> Matrix3x3:
    """Create a 3x3 identity matrix."""
    return np.identity(3)


def kronecker_delta(i: int, j: int) -> float:
    """
    Evaluates as 1 when parameters i and j are equal and 0 otherwise.

    :param i: First argument to check.
    :param j: Second argument to check.
    """
    return 1.0 if i == j else 0.0


def matrix_inverse(matrix: Matrix3x3Like) -> Matrix3x3:
    """
    Calculate the inverse :math:`M^{-1}` of a 3x3 matrix :math:`M`.

    The inverse :math:`M^{-1}` is the matrix which fufills the identity:

    .. math:: M^{-1} M = M M^{-1} = I

    where :math:`I` is the identity matrix.

    :param matrix: Matrix :math:`M^{-1}` to invert.
    :raises ValueError: Matrix is singular and hence an inverse does not exist.
    :return: Inverse of the matrix :math:`M`
    """
    try:
        return np.linalg.inv(matrix)  # type: ignore
    except np.linalg.LinAlgError as e:
        raise ValueError("Matrix is singular and hence cannot be inverted.") from e


def transpose(matrix: Matrix3x3Like) -> Matrix3x3:
    """Calculate the transpose of 3x3 matrices."""
    return np.swapaxes(matrix, -2, -1)
