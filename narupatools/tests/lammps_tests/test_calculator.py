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

import pytest

lammps = pytest.importorskip("lammps")

from narupatools.core.units import calorie, electronvolt, kilo, mole
from narupatools.lammps._converter import atoms_from_lammps_simulation
from narupatools.lammps._simulation import LAMMPSSimulation


@pytest.fixture(scope="module")
def simulation():
    return LAMMPSSimulation.from_file("./in.peptide")


@pytest.fixture
def atoms(simulation):
    return atoms_from_lammps_simulation(simulation)


def test_energy(atoms):
    # energy output by LAMMPS
    initial_energy = -6372.3759 * ((kilo * calorie / mole) >> (electronvolt))
    assert atoms.get_potential_energy() == pytest.approx(initial_energy, rel=1e-3)
