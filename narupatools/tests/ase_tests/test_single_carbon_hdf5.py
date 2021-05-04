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

import numpy as np
import pytest
from ase import Atoms
from ase.md import Langevin

from narupatools.ase import ASEDynamics, NullCalculator, UnitsASE
from narupatools.core import UnitsNarupa
from test_classes.single_carbon_hdf5 import SingleCarbonHDF5Tests

_NarupaToASE = UnitsNarupa >> UnitsASE


class TestASESingleCarbonHDF5(SingleCarbonHDF5Tests):
    @pytest.fixture
    def dynamics(self):
        atoms = Atoms(
            symbols=["C"],
            masses=np.array([12.000]) * _NarupaToASE.mass,
            positions=np.array([[5.0, 5.0, 5.0]]) * _NarupaToASE.length,
        )
        atoms.set_calculator(NullCalculator())
        langevin = Langevin(
            atoms,
            friction=0.01 / _NarupaToASE.time,
            fixcm=True,
            temperature_K=300,
            timestep=0.01 * _NarupaToASE.time,
        )
        return ASEDynamics(langevin)
