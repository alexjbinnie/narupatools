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

from narupatools.lammps.atom_properties import AtomID

lammps = pytest.importorskip("lammps")


from narupatools.lammps._simulation import LAMMPSSimulation


@pytest.fixture(scope="module")
def simulation():
    return LAMMPSSimulation.from_file("./in.peptide")


@pytest.fixture(scope="module")
def simulation_dynamic():
    return LAMMPSSimulation.from_file("./in.peptide")


def test_energy(simulation):
    assert simulation.potential_energy == pytest.approx(-26662.02362458588)


def test_temperature(simulation):
    assert simulation.temperature == pytest.approx(275.0, abs=15.0)


def test_positions(simulation):
    assert len(simulation.positions) == 2004
    assert simulation.positions[0] == pytest.approx(
        np.array([4.399993, 5.852678, 3.67855])
    )


def test_velocities(simulation):
    assert len(simulation.velocities) == 2004


def test_forces(simulation):
    assert len(simulation.forces) == 2004


def test_masses(simulation):
    assert len(simulation.masses) == 2004
    assert simulation.masses[0] > 0


def test_ids_sorted(simulation):
    assert np.all(np.diff(simulation.gather_atoms(AtomID)) >= 1)


def test_set_positions(simulation_dynamic):
    positions = simulation_dynamic.positions.copy()
    positions[0] += [10.0, 0.0, 0.0]
    simulation_dynamic.positions = positions
    new_positions = simulation_dynamic.positions
    assert new_positions[0] == pytest.approx(positions[0])
