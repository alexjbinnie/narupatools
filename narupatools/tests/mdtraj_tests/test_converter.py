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

import mdtraj
import pytest
from mdtraj import Trajectory
from narupa.trajectory import FrameData
from test_classes.converter import NeuraminidaseTestConverter

from narupatools.mdtraj.converter import (
    mdtraj_topology_to_frame,
    mdtraj_trajectory_to_frame,
)


@pytest.fixture(scope="module")
def mdtraj_trajectory(neuraminidase_pdb_filename) -> Trajectory:
    return mdtraj.load(neuraminidase_pdb_filename)


@pytest.fixture(scope="module")
def neuraminidase_frame(mdtraj_trajectory) -> FrameData:
    frame = FrameData()
    mdtraj_trajectory_to_frame(mdtraj_trajectory, frame=frame)
    mdtraj_topology_to_frame(mdtraj_trajectory.topology, frame=frame)
    return frame


@pytest.mark.converter
class TestMDTrajConverter(NeuraminidaseTestConverter):
    bond_count = 3321

    @pytest.fixture(autouse=True, scope="class")
    def frame(self, neuraminidase_frame):
        return neuraminidase_frame

    def test_import_works(self, mdtraj_trajectory):
        assert isinstance(mdtraj_trajectory, Trajectory)
