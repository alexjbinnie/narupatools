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

"""Methods for dealing with rotations and transformations."""

import numpy as np
from scipy.spatial.transform import Rotation

from narupatools.physics.rigidbody import center_of_mass
from narupatools.physics.typing import (
    Matrix3x3,
    Quaternion,
    QuaternionLike,
    RotationLike,
    ScalarArrayLike,
    Vector3Array,
    Vector3ArrayLike,
    Vector3Like,
)
from narupatools.physics.vector import sqr_magnitude


def quaternion_conjugate(quaternion: QuaternionLike, /) -> Quaternion:
    """
    Calculate the conjugate of a quaternion :math:`q = x i + y j + z k + w`.

    This is given by:

    .. math:: q^* = - x i - y j - z k + w

    :param quaternion: Quaternion :math:`q` to conjugate.
    :return: Conjugate :math:`q^*` of the quaternion :math:`q`.
    """
    return np.array([-quaternion[0], -quaternion[1], -quaternion[2], quaternion[3]])


def quaternion_inverse(quaternion: QuaternionLike, /) -> Quaternion:
    r"""
    Calculate the inverse of a quaternion :math:`q = x i + y j + z k + w`.

    This is given by:

    .. math:: q^{-1} = \frac{q^*}{|q|^2}

    where :math:`q^*` is the conjugate of the quaternion :math:`q` and `|q|` is the
    magnitude of the quaternion :math:`q`.

    :param quaternion: Quaternion :math:`q` to invert.
    :return: Inverse :math:`q^{-1}` of the quaternion :math:`q`.
    """
    return quaternion_conjugate(quaternion) / sqr_magnitude(quaternion)  # type: ignore


def quaternion_as_rotation_matrix(quaternion: QuaternionLike, /) -> Matrix3x3:
    r"""
    Calculate the rotation matrix from a quaternion :math:`q = x i + y j + z k + w`.

    This uses the formula:

    .. math:: R = \begin{bmatrix} 1 - 2(y^2 + z^2) & 2 (x y - z w) & 2 (x z + y w) \\
              2 (x y + z w) & 1 - 2(x^2 + z^2) & 2 (y z - x w) \\
              2(x z - y w) & 2 (y z + x w) & 1 - 2 (x^2 + y^2) \end{bmatrix}

    where :math:`q` is assumed to be a unit quaternion.

    :param quaternion: Quaternion :math:`q` to represent as a rotation matrix
    :return: Rotation matrix :math:`R` as a (3, 3) NumPy array.
    """
    return Rotation.from_quat(quaternion).as_matrix()  # type: ignore[arg-type]


def _as_rotation_matrix(rotation: RotationLike, /) -> Matrix3x3:
    rotation = np.asfarray(rotation)
    if rotation.shape == (3, 3):
        return rotation  # type: ignore
    elif rotation.shape == (4,):
        return quaternion_as_rotation_matrix(rotation)
    else:
        raise ValueError(f"Cannot interpret {rotation} as a rotation.")


def rotate_around_center_of_mass(
    *, masses: ScalarArrayLike, positions: Vector3ArrayLike, rotation: RotationLike
) -> Vector3Array:
    """
    Rotate a set of particles about their center of mass.

    :param masses: List of masses :math:`m_i` of each particle.
    :param positions: List of positions :math:`r_i` of each particle.
    :param rotation: Rotation to perform, either as a 3x3 rotation matrix :math:`R` or
                     a quaternion :math:`q = x i + y j + z k + w` as an array
                     (x, y, z, w).
    :return: Set of positions that have been rotated around their center of mass.
    """
    return rotate_around_point(
        positions=positions,
        rotation=rotation,
        origin=center_of_mass(masses=masses, positions=positions),
    )


def rotate_around_point(
    *, positions: Vector3ArrayLike, rotation: RotationLike, origin: Vector3Like
) -> Vector3Array:
    """
    Rotate a set of points about an arbitrary point :math:`c`.

    :param origin: Point :math:`c` about which to rotate the positions.
    :param positions: List of positions :math:`r_i`.
    :param rotation: Rotation to perform, either as a 3x3 rotation matrix :math:`R` or
                     a quaternion :math:`q = x i + y j + z k + w` as an array
                     (x, y, z, w).
    :return: Set of positions that have been rotated around the origin :math:`c`.
    """
    rotation = _as_rotation_matrix(rotation)
    origin = np.asfarray(origin)
    return np.array(
        list(map(lambda p: origin + np.matmul(rotation, p - origin), positions))
    )


def rotate_around_origin(
    *, positions: Vector3ArrayLike, rotation: RotationLike
) -> Vector3Array:
    """
    Rotate a set of points about the origin (0, 0, 0).

    :param positions: List of positions :math:`r_i`.
    :param rotation: Rotation to perform, either as a 3x3 rotation matrix :math:`R` or
                     a quaternion :math:`q = x i + y j + z k + w` as an array
                     (x, y, z, w).
    :return: Set of positions that have been rotated around the origin (0, 0, 0).
    """
    rotation = _as_rotation_matrix(rotation)
    return np.array(list(map(lambda p: np.matmul(rotation, p), positions)))
