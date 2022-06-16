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

"""Utility methods for using vectors."""

from typing import Optional, Union, overload

import numpy as np
import numpy.typing as npt

from ._quaternion import quaternion
from .typing import Matrix3x3, Vector3, Vector3Like, VectorN, VectorNLike


def vector(*args: Union[float, int]) -> VectorN:
    """
    Create a vector from a set of coordinates.

    This creates a NumPy array from the arguments with its dtype set to float. This is
    intended as a shorthand for creating vectors which are guaranteed to be floats and
    not integers.

    :param args: Components of the vector.
    """
    return np.array(args, dtype=float)


def dot_product(a: VectorN, b: VectorN, /) -> float:
    r"""
    Calculate the dot product of two N-dimensional vectors :math:`a` and :math:`b`.

    This is defined as:

    .. math:: a \cdot b = \sum_i a_i b_i

    If both a and b are one dimensional vectors, the dot product is calculated as above.

    If either a or b have more dimensions, NumPy broadcasting is used:

    * (N) + (N) -> float
    * (M, N) + (N) -> (M)
    * (M, N) + (M, N) -> (M)
    * (M, N) + (K, M, N) -> (K, M)

    :param a: Vector :math:`a`.
    :param b: Vector :math:`b`.
    :return: Dot product :math:`a \cdot b` of the vectors :math:`a` and :math:`b`.
    """
    a = np.asfarray(a)
    b = np.asfarray(b)
    return (a * b).sum(axis=-1)  # type: ignore


def cross_product(a: Vector3Like, b: Vector3Like, /) -> Vector3:
    r"""
    Calculate the cross product of two 3-dimensional vectors :math:`a` and :math:`b`.

    This is defined as:

    .. math:: a \times b = \begin{vmatrix} \hat i & \hat j & \hat k \\ a_x & a_y & a_z
              \\ b_x & b_y & b_z \end{vmatrix}

    :param a: Vector :math:`a`.
    :param b: Vector :math:`b`.
    :return: Cross product :math:`a \times b` of the vectors :math:`a` and :math:`b`.
    """
    a = np.asfarray(a)
    b = np.asfarray(b)
    if a.shape == (3,) and b.shape == (3,):
        return np.array(
            [
                a[1] * b[2] - b[1] * a[2],
                a[2] * b[0] - b[2] * a[0],
                a[0] * b[1] - b[0] * a[1],
            ]
        )
    else:
        return np.cross(a, b)


def zero_vector() -> Vector3:
    """Zero vector in three dimensions."""
    return np.zeros(3)


def sqr_magnitude(vector: VectorNLike, /) -> float:
    """Get the square magnitude of a n-dimensional vector."""
    return (vector * vector).sum(axis=-1)  # type: ignore


def magnitude(vector: VectorNLike, /) -> float:
    """Get the magnitude of a n-dimensional vector."""
    return np.linalg.norm(vector, axis=-1)  # type: ignore


@overload
def normalized(vector: quaternion, /) -> quaternion:
    ...


@overload
def normalized(vector: npt.NDArray[quaternion], /) -> npt.NDArray[quaternion]:
    ...


@overload
def normalized(vector: VectorNLike, /) -> npt.NDArray[np.floating]:
    ...


def normalized(
    vector: Union[quaternion, npt.NDArray[quaternion], VectorNLike], /
) -> Union[quaternion, npt.NDArray[quaternion], npt.NDArray[np.floating]]:
    """Normalize an n-dimensional vector."""
    if isinstance(vector, quaternion):
        return np.normalized(vector)  # type: ignore
    arr = np.asarray(vector)
    if arr.dtype == quaternion:
        return np.normalized(arr)  # type: ignore
    else:
        mag = np.sqrt((arr**2).sum(axis=-1))
        if isinstance(mag, float):
            return np.nan_to_num(arr / mag)
        else:
            return np.nan_to_num(arr / mag[:, np.newaxis])  # type: ignore


def vector_projection(vector: Vector3Like, onto: Vector3Like, /) -> Vector3:
    r"""
    Calculate the **vector projection** :math:`a_1` of :math:`a` onto :math:`b`.

    This is given by:

    .. math:: a_1 = \frac{a \cdot b}{b \cdot b} b

    The vector projection gives a vector in the same direction as b.

    This is also known as the vector component or vector resolution of :math:`a` in the
    direction of :math:`b`.

    :param vector: Vector :math:`a` to project.
    :param onto: Vector :math:`b` to project onto.
    :return: Vector projection :math:`a_1` of vector :math:`a` onto vector :math:`b`.
    """
    vector = np.asfarray(vector)
    onto = np.asfarray(onto)
    dot = (onto * onto).sum(axis=-1)
    return ((vector * onto).sum(axis=-1) / dot)[..., np.newaxis] * onto  # type: ignore


def vector_rejection(vector: Vector3Like, onto: Vector3Like, /) -> Vector3:
    r"""
    Calculate the **vector rejection** :math:`a_2` of :math:`a` from :math`b`.

    This is given by:

    .. math:: a_2 = a - a_1 = a - \frac{a \cdot b}{b \cdot b} b

    where :math:`a_1` is the vector projection of :math:`a` onto :math:`b`.

    The vector rejection is the component of :math:`a` which is perpendicular to
    :math:`b`.

    :param vector: Vector :math:`a` to reject.
    :param onto: Vector :math:`b` to reject from.
    :return: Vector rejection :math:`a_2` of vector :math:`a` from vector :math:`b`.
    """
    return np.asfarray(vector) - vector_projection(vector, onto)


def distance(vector1: Vector3Like, vector2: Optional[Vector3Like] = None, /) -> float:
    r"""
    Calculate the distance :math:`d` between two points :math:`a` and :math:`b`.

    :param vector1: Point :math:`a`.
    :param vector2: Point :math:`b`.
    :return: Distance between the two points.
    """
    vector1 = np.asfarray(vector1)
    if vector2 is None:
        if vector1.shape[-2] != 2:
            raise ValueError(
                "Cannot take distances of array without second last axes being shape 2."
            )
        return np.linalg.norm(vector1[..., 1, :] - vector1[..., 0, :], axis=-1)  # type: ignore
    return np.sqrt(((vector1[..., :] - vector2[..., :]) ** 2).sum(axis=-1))  # type: ignore


def sqr_distance(
    point1: Vector3Like, point2: Vector3Like, /
) -> Union[float, np.ndarray]:
    r"""
    Calculate the square distance between two points :math:`a` and :math:`b`.

    :param point1: Point :math:`a`.
    :param point2: Point :math:`b`.
    :return: Square distance between the two points.
    """
    return ((point1[..., :] - point2[..., :]) ** 2).sum(axis=-1)  # type: ignore


def angle(vector1: Vector3Like, vector2: Vector3Like, /) -> Union[float, np.ndarray]:
    r"""
    Calculate the angle in radians between two vectors :math:`a` and  :math:`b`.

    This uses the relation:

    .. math:: a \cdot b = |a| |b| \cos \theta

    :param vector1: Vector :math:`a`.
    :param vector2: Vector :math:`b`.
    :raises ValueError: One of the vectors is zero.
    :return: Angle between the two vectors in radians.
    """
    a = np.asfarray(vector1)
    b = np.asfarray(vector2)
    return np.arccos(  # type: ignore
        (a * b).sum(axis=-1) / np.sqrt((a * a).sum(axis=-1) * (b * b).sum(axis=-1))
    )


def cross_product_matrix(vector: Vector3Like, /) -> Matrix3x3:
    r"""
    Calculate the matrix `F` that represents the left cross product with :math:`a`.

    This matrix fulfills for all vectors :math:`v` the identity:

    .. math:: F v = a \times v

    This converts taking the cross product on the left by :math:`a` to a matrix
    multiplication.

    :param vector: Vector :math:`a`.
    :return: Skew-symmetric matrix :math:`F` that satisfies :math:`F v = a \times v` for
             all :math:`v`.
    """
    return np.array(
        [
            [0, -vector[2], vector[1]],
            [vector[2], 0, -vector[0]],
            [-vector[1], vector[0], 0],
        ],
        dtype=float,
    )


def right_cross_product_matrix(vector: Vector3Like, /) -> Matrix3x3:
    r"""
    Calculate the matrix that represents the right cross product with :math:`a`.

    This matrix fulfills for all vectors :math:`v` the identity:

    .. math:: F v = v \times a

    This converts taking the cross product on the right by :math:`a` to a matrix
    multiplication.

    :param vector: Vector :math:`a`.
    :return: Skew-symmetric matrix :math:`F` that satisfies :math:`F v = v \times a` for
             all :math:`v`.
    """
    return np.array(
        [
            [0, vector[2], -vector[1]],
            [-vector[2], 0, vector[0]],
            [vector[1], -vector[0], 0],
        ],
        dtype=float,
    )


def left_vector_triple_product_matrix(a: Vector3Like, b: Vector3Like, /) -> Matrix3x3:
    r"""
    Calculate the matrix form of the left cross product with :math:`a` and :math:`b`.

    This matrix fulfills for all vectors :math:`v` the identity:

    .. math:: F v = a \times (b \times v)

    This converts taking the cross product on the left by :math:`b` and then by
    :math:`a` to a matrix multiplication.

    :param a: Vector :math:`a`.
    :param b: Vector :math:`b`.
    """
    return np.matmul(cross_product_matrix(a), cross_product_matrix(b))  # type: ignore


def outer_product(a: Vector3Like, b: Vector3Like, /) -> np.ndarray:
    r"""
    Calculate the outer product of two vectors :math:`a` and :math:`b`.

    This is a matrix defined as:

    .. math:: A = \begin{pmatrix} a_1 b_1 & \cdots & a_1 b_3 \\ \vdots & \ddots &
              \vdots \\ a_3 b_1 & \cdots & a_3 b_3 \end{pmatrix}

    :param a: Vector :math:`a`.
    :param b: Vector :math:`b`.
    :return: Matrix formed by the outer product of :math:`a` and :math:`b`.
    """
    return np.einsum("...i,...j->...ij", a, b)  # type: ignore
