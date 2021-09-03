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

"""Methods for generating random values."""

import math
import random

import numpy as np

from narupatools.core.random import random_float, random_integer
from narupatools.physics.typing import Vector3

from .quaternion import quaternion


def random_sphere() -> Vector3:
    """Generate a random point on the unit sphere."""
    theta = 2 * math.pi * random.random()
    phi = math.acos(2 * random.random() - 1)
    return np.array(
        [
            math.cos(theta) * math.sin(phi),
            math.sin(theta) * math.sin(phi),
            math.cos(phi),
        ]
    )


def random_vector(*, max_magnitude: float = 1.0) -> Vector3:
    """
    Generate a random vector of length 3.

    :param max_magnitude: Maximum magnitude of the vector.
    """
    return random_sphere() * max_magnitude


def random_quaternion() -> quaternion:
    """Generate a random vector of length 3."""
    return quaternion(*[random_float(min=-1.0, max=1.0) for _ in range(4)])


def random_unit_quaternion() -> quaternion:
    """Generate a random vector of length 3."""
    return np.normalized(random_quaternion())
