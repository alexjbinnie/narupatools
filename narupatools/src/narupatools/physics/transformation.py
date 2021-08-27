"""Spatial transformations."""

from __future__ import annotations

from typing import ClassVar, Union, overload

import numpy as np

from narupatools.physics.rigidbody import center_of_mass
from narupatools.physics.typing import (
    Matrix3x3,
    ScalarArrayLike,
    Vector3,
    Vector3Array,
    Vector3ArrayLike,
    Vector3Like,
)
from narupatools.physics.vector import vector

from .quaternion import (
    as_rotation_matrix,
    as_rotation_vector,
    from_rotation_vector,
    quaternion,
)


class Translation:
    """
    Translation in three dimensional space.

    A translation shifts all points in space by a three dimensional vector :math:`v`.
    """

    identity: ClassVar[Translation]
    """Identity transformation."""

    def __init__(self, translation: np.ndarray):
        self.__vector = translation

    @property
    def translation(self) -> Vector3:
        """Translation vector of this transformation."""
        return self.__vector

    @overload
    def __mul__(self, other: Translation) -> Translation:
        ...

    @overload
    def __mul__(self, other: np.ndarray) -> np.ndarray:
        ...

    def __mul__(
        self, other: Union[Translation, np.ndarray]
    ) -> Union[Translation, np.ndarray]:
        if isinstance(other, Translation):
            return Translation(self.__vector + other.__vector)
        if (
            isinstance(other, np.ndarray)
            and len(other.shape) == 2
            and other.shape[1] == 3
        ):
            return np.array([self.__vector + vector for vector in other])
        return NotImplemented

    def __add__(self, other: Translation) -> Translation:
        if isinstance(other, Translation):
            return Translation(self.__vector + other.__vector)
        return NotImplemented  # type: ignore[unreachable]

    def __invert__(self) -> Translation:
        return Translation(-self.__vector)

    def __array__(self) -> np.ndarray:
        return self.__vector


Translation.identity = Translation(vector(0, 0, 0))


class Rotation:
    """
    Represents a rotation in 3-dimensional space.

    Internally, rotations are stored as quaternions, due to their advantages over other forms such as
    rotation matrices, rotation vectors and euler angles.
    """

    identity: ClassVar[Rotation]
    """Identity transformation."""

    def __init__(self, quat: quaternion):
        self.__quat = quat
        """Quaternion stored in the form (qx i + qy j + qz k + qw), with the scalar part last."""

    @property
    def versor(self) -> quaternion:
        """Unit quaternion (versor) representing this rotation."""
        return self.__quat

    @property
    def inverse(self) -> Rotation:
        """Inverse of this rotation."""
        return ~self

    @property
    def rotation_vector(self) -> Vector3:
        """Rotation vector representing this rotation."""
        return as_rotation_vector(self.__quat)

    @property
    def rotation_matrix(self) -> Matrix3x3:
        """Rotation matrix representing this rotation."""
        return as_rotation_matrix(self.__quat)

    def __invert__(self) -> Rotation:
        return Rotation(self.__quat.conjugate())

    @overload
    def __mul__(self, other: Rotation) -> Rotation:
        ...

    @overload
    def __mul__(self, other: np.ndarray) -> np.ndarray:
        ...

    def __mul__(
        self, other: Union[Rotation, np.ndarray]
    ) -> Union[Rotation, np.ndarray]:
        if isinstance(other, Rotation):
            return Rotation(self.__quat * other.__quat)
        if other.shape == (3,):
            return np.matmul(self.rotation_matrix, other)  # type: ignore[no-any-return]
        elif len(other.shape) == 2 and other.shape[1] == 3:
            return np.dot(self.rotation_matrix, other.T).T  # type: ignore[no-any-return]
        return NotImplemented

    @classmethod
    def from_rotation_vector(cls, rot: Vector3) -> Rotation:
        """Create a rotation from a rotation vector."""
        rot = np.array(rot, dtype=float)
        if rot.shape != (3,):
            raise ValueError(f"Rotation vector has invalid shape {rot.shape}")
        return Rotation(from_rotation_vector(rot))

    def __array__(self) -> np.ndarray:
        return self.__quat.components

    def rotate_around_center_of_mass(
        self,
        *,
        masses: ScalarArrayLike,
        positions: Vector3ArrayLike,
    ) -> Vector3Array:
        """
        Rotate a set of particles about their center of mass.

        :param masses: List of masses :math:`m_i` of each particle.
        :param positions: List of positions :math:`r_i` of each particle.
        :return: Set of positions that have been rotated around their center of mass.
        """
        return self.rotate_around_point(
            positions=positions,
            origin=center_of_mass(masses=masses, positions=positions),
        )

    def rotate_around_point(
        self, *, positions: Vector3ArrayLike, origin: Vector3Like
    ) -> Vector3Array:
        """
        Rotate a set of points about an arbitrary point :math:`c`.

        :param origin: Point :math:`c` about which to rotate the positions.
        :param positions: List of positions :math:`r_i`.
        :return: Set of positions that have been rotated around the origin :math:`c`.
        """
        rotation_matrix = self.rotation_matrix
        origin = np.asfarray(origin)
        return np.array(
            list(
                map(
                    lambda p: origin + np.matmul(rotation_matrix, p - origin), positions
                )
            )
        )


Rotation.identity = Rotation(quaternion(1, 0, 0, 0))
