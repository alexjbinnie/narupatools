from ase.md import Langevin

from narupatools.ase.openmm import ASEOpenMMDynamics


def test_from_path(nanotube_xml_filename):
    dynamics = ASEOpenMMDynamics.from_xml_file(nanotube_xml_filename)
    assert len(dynamics.atoms) == 65
    assert isinstance(dynamics.molecular_dynamics, VelocityVerlet)


def test_from_string(nanotube_xml_filename):
    with open(nanotube_xml_filename) as file:
        contents = file.read()
    dynamics = ASEOpenMMDynamics.from_xml_string(contents)
    assert len(dynamics.atoms) == 65
    assert isinstance(dynamics.molecular_dynamics, VelocityVerlet)
