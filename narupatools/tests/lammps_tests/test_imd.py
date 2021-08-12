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

from narupatools.lammps import LAMMPSSimulation


@pytest.fixture
def simulation():
    simulation = LAMMPSSimulation.from_file("./in.peptide")
    simulation.add_imd_force()
    return simulation


def test_initial_no_imd(simulation):
    assert np.all(simulation.get_imd_forces() == 0)


def test_set_imd_forces(simulation):
    simulation.set_imd_force(2, np.array([10.0, 0.0, 0.0]))
    simulation.run(0)
    assert simulation.get_imd_forces()[2] == pytest.approx([10.0, 0.0, 0.0])


def test_imd_forces_changes_actual_force(simulation):
    existing_force = simulation.forces[2]
    simulation.set_imd_force(2, np.array([100.0, 0.0, 0.0]))
    simulation.run(0)
    assert simulation.forces[2] - existing_force == pytest.approx(
        [100.0, 0.0, 0.0], rel=1e-3, abs=0.1
    )
