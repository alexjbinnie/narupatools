import pytest
from openmm import LangevinIntegrator, VerletIntegrator

from narupatools.openmm import OpenMMSimulation, deserialize_simulation


@pytest.fixture
def nanotube_openmm_simulation(nanotube_xml_filename):
    with open(nanotube_xml_filename) as file:
        return deserialize_simulation(file.read())


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
