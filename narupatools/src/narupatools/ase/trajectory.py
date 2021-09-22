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
#
# Originally part of the narupa-ase package.
# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Modified under the terms of the GPL.

"""Implements a set of ASE atoms objects as a trajectory that can be played back."""

from typing import Optional, Sequence

from ase.atoms import Atoms
from infinite_sets import InfiniteSet
from MDAnalysis import Universe
from narupa.trajectory import FrameData

from narupatools.ase import ase_atoms_to_frame
from narupatools.frame.frame_source import TrajectorySource
from narupatools.mdanalysis import mdanalysis_universe_to_frame


class ASETrajectory(TrajectorySource):
    """Trajectory playback using one ore more  ASE `Atoms` objects."""

    def __init__(
        self, trajectory: Sequence[Atoms], *, universe: Optional[Universe] = None
    ):
        self._trajectory = trajectory
        self._universe = universe

    def __len__(self) -> int:
        return len(self._trajectory)

    def get_frame(  # noqa: D102
        self, *, index: int, fields: InfiniteSet[str]
    ) -> FrameData:
        frame = FrameData()
        if self._universe:
            ase_atoms_to_frame(self._trajectory[index], fields=fields, frame=frame)
            added_fields = set(frame.arrays.keys()) | set(frame.values.keys())
            mdanalysis_universe_to_frame(
                self._universe, fields=fields - added_fields, frame=frame
            )
            return frame
        else:
            ase_atoms_to_frame(self._trajectory[index], fields=fields, frame=frame)
        return frame
