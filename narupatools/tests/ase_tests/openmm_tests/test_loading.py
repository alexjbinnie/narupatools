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

from ase.md import VelocityVerlet

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
