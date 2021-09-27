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
from simtk.openmm import LangevinIntegrator
from simtk.openmm.app import Simulation
from simtk.unit import kelvin, picosecond, picoseconds

from narupatools.ase._units import UnitsASE
from narupatools.ase.openmm._converter import (
    openmm_simulation_to_ase_molecular_dynamics,
)
from narupatools.physics.units import UnitsNarupa


@pytest.fixture(params=[0.01, 0.1, 1.0])
def timestep(request):
    return request.param


@pytest.fixture(params=[50, 300, 1000])
def temperature(request):
    return request.param


@pytest.fixture(params=[1e-2, 1e-1, 1.0])
def friction(request):
    return request.param


@pytest.fixture
def make_villin_openmm_langevin_simulation(villin_openmm_pdbfile, villin_openmm_system):
    def make(temperature, friction, timestep):
        integrator = LangevinIntegrator(
            temperature * kelvin, friction / picosecond, timestep * picoseconds
        )
        simulation = Simulation(
            villin_openmm_pdbfile.topology, villin_openmm_system, integrator
        )
        simulation.context.setPositions(villin_openmm_pdbfile.positions)
        return simulation

    return make


@pytest.fixture
def villin_ase_langevin(villin_openmm_langevin):
    return openmm_simulation_to_ase_molecular_dynamics(villin_openmm_langevin)


ASEToNarupa = UnitsASE >> UnitsNarupa


class TestOpenMMToASEConverter:
    def test_temperature_read_from_openmm(
        self, make_villin_openmm_langevin_simulation, temperature
    ):
        dynamics = make_villin_openmm_langevin_simulation(temperature, 1.0, 1.0)
        ase_dynamics = openmm_simulation_to_ase_molecular_dynamics(dynamics)
        assert ase_dynamics.todict()["temperature_K"] == pytest.approx(temperature)

    def test_friction_read_from_openmm(
        self, make_villin_openmm_langevin_simulation, friction
    ):
        dynamics = make_villin_openmm_langevin_simulation(300, friction, 1.0)
        ase_dynamics = openmm_simulation_to_ase_molecular_dynamics(dynamics)
        assert ase_dynamics.todict()["friction"] / ASEToNarupa.time == pytest.approx(
            friction
        )

    def test_timestep_read_from_openmm(
        self, make_villin_openmm_langevin_simulation, timestep
    ):
        dynamics = make_villin_openmm_langevin_simulation(300, 1.0, timestep)
        ase_dynamics = openmm_simulation_to_ase_molecular_dynamics(dynamics)
        assert ase_dynamics.todict()["timestep"] * ASEToNarupa.time == pytest.approx(
            timestep
        )
