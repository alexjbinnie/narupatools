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
from openmm import LangevinIntegrator, VerletIntegrator

from narupatools.openmm import OpenMMSimulation, deserialize_simulation


@pytest.fixture
def nanotube_openmm_simulation(nanotube_xml_filename):
    with open(nanotube_xml_filename) as file:
        return deserialize_simulation(file.read())


def test_platform(nanotube_openmm_simulation):
    print(nanotube_openmm_simulation.context.getPlatform().getName())


@pytest.fixture
def nanotube_simulation(nanotube_openmm_simulation):
    return OpenMMSimulation.from_simulation(nanotube_openmm_simulation)


def test_simulation_not_copy(nanotube_openmm_simulation):
    simulation = OpenMMSimulation.from_simulation(nanotube_openmm_simulation)

    assert simulation.system is nanotube_openmm_simulation.system
    assert simulation.topology is nanotube_openmm_simulation.topology
    assert simulation.integrator is nanotube_openmm_simulation.integrator
    assert simulation.context is nanotube_openmm_simulation.context


def test_change_integrator(nanotube_simulation):
    assert isinstance(nanotube_simulation.context.getIntegrator(), VerletIntegrator)

    nanotube_simulation.integrator = LangevinIntegrator(300, 10, 0.001)

    assert isinstance(nanotube_simulation.context.getIntegrator(), LangevinIntegrator)


def test_global_parameters_len(nanotube_simulation):
    assert len(nanotube_simulation.global_parameters) == 6


def test_global_parameters_keys(nanotube_simulation):
    assert nanotube_simulation.global_parameters.keys() == {
        "c12",
        "c6",
        "cs",
        "mx",
        "my",
        "mz",
    }


def test_global_parameters_values(nanotube_simulation):
    assert list(nanotube_simulation.global_parameters.values()) == [
        1.0,
        0.0,
        -25.5,
        7.0,
        7.0,
        7.0,
    ]


def test_global_parameters_contains(nanotube_simulation):
    assert "my" in nanotube_simulation.global_parameters


def test_global_parameters_not_contains(nanotube_simulation):
    assert "not_key" not in nanotube_simulation.global_parameters


def test_global_parameters_set(nanotube_simulation):
    nanotube_simulation.global_parameters["c6"] = 2.0


def test_global_parameters_set_not_defined(nanotube_simulation):
    with pytest.raises(KeyError):
        nanotube_simulation.global_parameters["not_key"] = 0.0
