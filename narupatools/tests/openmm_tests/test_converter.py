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
from narupa.trajectory import FrameData
from simtk.openmm.app import PDBFile, Topology

from narupatools.frame import ParticlePositions
from narupatools.openmm import openmm_topology_to_frame
from test_classes.converter import NeuraminidaseTestConverter


@pytest.fixture(scope="module")
def neuraminidase_pdbfile(neuraminidase_pdb_filename) -> PDBFile:
    return PDBFile(neuraminidase_pdb_filename)


@pytest.fixture(scope="module")
def neuraminidase_frame(neuraminidase_pdbfile) -> FrameData:
    frame = FrameData()
    openmm_topology_to_frame(neuraminidase_pdbfile.getTopology(), existing=frame)
    ParticlePositions.set(frame, neuraminidase_pdbfile.getPositions(asNumpy=True))
    return frame


@pytest.mark.converter
class TestOpenMMConverter(NeuraminidaseTestConverter):
    @pytest.fixture(autouse=True, scope="class")
    def frame(self, neuraminidase_frame):
        return neuraminidase_frame

    def test_import_works(self, neuraminidase_pdbfile):
        assert isinstance(neuraminidase_pdbfile, PDBFile)
        assert isinstance(neuraminidase_pdbfile.getTopology(), Topology)
