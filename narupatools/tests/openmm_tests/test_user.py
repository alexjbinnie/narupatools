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

from narupatools.ase.openmm import ASEOpenMMDynamics
from narupatools.openmm import OpenMMDynamics


def test_ase_openmm(neuraminidase_openmm_xml_filename):
    dynamics = ASEOpenMMDynamics.from_xml_file(neuraminidase_openmm_xml_filename)
    dynamics.run(100)


def test_openmm(neuraminidase_openmm_xml_filename):
    dynamics = OpenMMDynamics.from_xml_file(neuraminidase_openmm_xml_filename)
    dynamics.run(100)
