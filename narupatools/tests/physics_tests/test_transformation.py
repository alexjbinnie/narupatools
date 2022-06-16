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
import math

import numpy as np
import pytest

from narupatools.physics.random import random_quaternion, random_vector
from narupatools.physics.transformation import Rotation
from narupatools.physics.vector import vector


@pytest.fixture
def rotation_vector(seed):
    return random_vector(max_magnitude=3.14)


@pytest.fixture
def unit_quaternion(seed):
    return random_quaternion().normalized()


def test_rotation_invertible(unit_quaternion):
    rotationa = Rotation(unit_quaternion)
    rotationb = ~rotationa
    rotationc = rotationa.inverse
    assert rotationb == pytest.approx(rotationc)
    assert (rotationa @ rotationb) == pytest.approx(Rotation.identity)
    assert (rotationb @ rotationa) == pytest.approx(Rotation.identity)


def test_rotation_vector_inverse(rotation_vector):
    rotation = Rotation.from_rotation_vector(rotation_vector)
    assert rotation.rotation_vector == pytest.approx(rotation_vector)
    reverse = (~rotation).rotation_vector
    assert reverse == pytest.approx(-rotation_vector)


def test_rotate_single_vector():
    point = vector(0, 1, 0)
    rot_vec = vector(0.5 * math.pi, 0.0, 0.0)
    rotation = Rotation.from_rotation_vector(rot_vec)
    assert rotation @ point == pytest.approx(vector(0, 0, 1))


def test_rotate_vector_array():
    points = np.array([vector(0, 1, 0), vector(0, 0, 1)])
    rot_vec = vector(0.5 * math.pi, 0.0, 0.0)
    rotation = Rotation.from_rotation_vector(rot_vec)
    assert rotation @ points == pytest.approx(
        np.array([vector(0, 0, 1), vector(0, -1, 0)])
    )
