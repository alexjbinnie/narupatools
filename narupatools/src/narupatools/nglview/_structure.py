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

"""Classes derived from Structure that allow NGLView to render new systems."""

from io import StringIO

import nglview
from ase.atoms import Atoms
from ase.io import write
from infinite_sets import everything
from narupa.trajectory.frame_data import FrameData

from narupatools.frame._converter import frame_to_pdb_string
from narupatools.frame._frame_source import FrameSource, TrajectorySource
from narupatools.frame.fields import ParticlePositions
from narupatools.physics.typing import Vector3Array


class ASEStructure(nglview.Structure):
    """ASE Structure for nglview that does not use temporary files."""

    def __init__(self, atoms: Atoms, /) -> None:
        super().__init__()
        self.atoms = atoms

    def get_structure_string(self) -> str:
        """Create a PDB string so NGLView can read in the structure."""
        file = StringIO("")
        write(file, self.atoms, "proteindatabank")
        return file.getvalue()


class FrameDataStructure(nglview.Structure):
    """FrameData Structure for nglview that does not use temporary files."""

    def __init__(self, frame: FrameData, /) -> None:
        super().__init__()
        self.frame = frame

    def get_structure_string(self) -> str:
        """Create a PDB string so NGLView can read in the structure."""
        return frame_to_pdb_string(self.frame)


class NarupaToolsFrame(nglview.Structure):
    """Wrapper around a FrameSource for use with nglview."""

    def __init__(self, frame: FrameSource):
        super().__init__()
        self.frame = frame

    def get_structure_string(self) -> str:  # noqa: D102
        frame = self.frame.get_frame(fields=everything())
        return frame_to_pdb_string(frame)


class NarupaToolsTrajectory(nglview.Trajectory, nglview.Structure):
    """Wrapper around a TrajectorySource for use with nglview."""

    def __init__(self, trajectory: TrajectorySource):
        super(nglview.Structure, self).__init__()
        super(nglview.Trajectory, self).__init__()
        self.trajectory = trajectory

    def get_coordinates(self, index: int) -> Vector3Array:  # noqa: D102
        frame = self.trajectory.get_frame(index=index, fields={ParticlePositions.key})
        positions = ParticlePositions.get(frame) * 10.0
        return positions

    @property
    def n_frames(self) -> int:  # noqa: D102
        return len(self.trajectory)

    def get_structure_string(self) -> str:  # noqa: D102
        frame = self.trajectory.get_frame(index=0, fields=everything())
        return frame_to_pdb_string(frame)
