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

from narupatools.physics.constants import boltzmann_constant
from narupatools.physics.thermodynamics import maxwell_boltzmann_velocities


def test_maxwell_boltzmann_mean():
    mass = 2.3
    temperature = 320

    masses = np.array([mass] * 1000000)

    velocities = maxwell_boltzmann_velocities(masses=masses, temperature=temperature)

    assert np.linalg.norm(velocities, axis=1).mean() == pytest.approx(
        np.sqrt(8 * boltzmann_constant * temperature / (math.pi * mass)), rel=1e-2
    )


def test_maxwell_boltzmann_kinetic_energy():
    masses = 0.1 + 10.0 * np.random.random_sample(1000000)
    temperature = 320

    ke = (
        0.5
        * masses
        * np.linalg.norm(
            maxwell_boltzmann_velocities(masses=masses, temperature=temperature), axis=1
        )
        ** 2
    ).mean()

    assert ke == pytest.approx(1.5 * boltzmann_constant * temperature, rel=1e-2)
