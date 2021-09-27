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

"""Methods relating to thermodynamics."""

import numpy as np
from numpy.random import standard_normal

from .constants import boltzmann_constant
from .typing import ScalarArray, Vector3Array


def maxwell_boltzmann_velocities(
    masses: ScalarArray, temperature: float
) -> Vector3Array:
    """
    Create a set of velocities that are distributed according to the Maxwell-Boltzmann distribution.

    :param masses: Masses of each particle in daltons.
    :param temperature: Thermodynamic temperature in kelvin.
    :return: NumPy array of velocities in nanometers per picoseconds.
    """
    dirs = standard_normal((len(masses), 3))
    return dirs * np.sqrt(boltzmann_constant * temperature / masses)[:, np.newaxis]  # type: ignore
