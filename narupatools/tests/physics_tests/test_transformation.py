import pytest

from narupatools.physics.matrix import matrix_inverse
from narupatools.physics.random import random_quaternion
from narupatools.physics.rotations import (
    quaternion_as_rotation_matrix,
    quaternion_inverse,
)


@pytest.fixture
def quaternion(seed):
    return random_quaternion()


def test_rotation_matrix_orthogonal(quaternion):
    rotation_matrix = quaternion_as_rotation_matrix(quaternion)
    assert matrix_inverse(rotation_matrix) == pytest.approx(rotation_matrix.T)


def test_rotation_matrix_quaternion_inverse(quaternion):
    rotation_matrix = quaternion_as_rotation_matrix(quaternion)
    rotation_matrix_inverse = quaternion_as_rotation_matrix(
        quaternion_inverse(quaternion)
    )
    assert matrix_inverse(rotation_matrix) == pytest.approx(rotation_matrix_inverse)
