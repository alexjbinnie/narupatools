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

from test_classes.single_carbon_hdf5 import SingleCarbonHDF5Tests

from narupatools.ase import ASEDynamics, UnitsASE
from narupatools.core import UnitsNarupa
from narupatools.lammps import (
    LAMMPSDynamics,
    LAMMPSSimulation,
    atoms_from_lammps_simulation,
)
from narupatools.lammps.regions import Box
from narupatools.physics.vector import vector

_NarupaToASE = UnitsNarupa >> UnitsASE


@pytest.fixture
def single_carbon_simulation():
    simulation = LAMMPSSimulation.create_new("real")
    simulation.command("atom_style atomic")
    simulation.command("boundary s s s")
    simulation.create_box(1, Box.bounds(vector(0, 0, 0), vector(10, 10, 10)))
    simulation.set_mass(type=1, mass=12.000)
    simulation.create_atom(type=1, position=vector(5, 5, 5))
    simulation.timestep = 0.01
    simulation.run(0)
    return simulation


@pytest.fixture
def single_carbon_langevin_simulation(single_carbon_simulation):
    single_carbon_simulation.command("fix 1 all nve")
    single_carbon_simulation.setup_langevin(temperature=1, friction=0.01, seed=22)
    single_carbon_simulation.run(0)
    return single_carbon_simulation


class TestLAMMPSSingleCarbonHDF5(SingleCarbonHDF5Tests):
    @pytest.fixture
    def dynamics(self, single_carbon_langevin_simulation):
        return LAMMPSDynamics(single_carbon_langevin_simulation)


class TestASELAMMPSSingleCarbonHF5(SingleCarbonHDF5Tests):
    @pytest.fixture
    def dynamics(self, single_carbon_simulation):
        atoms = atoms_from_lammps_simulation(single_carbon_simulation)
        dynamics = ASEDynamics.create_langevin(
            atoms, friction=0.01, timestep=0.01, temperature=300
        )
        return dynamics
