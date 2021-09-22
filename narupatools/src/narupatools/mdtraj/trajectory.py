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

"""Adapts a MDTraj trajectory for trajectory playback."""

from __future__ import annotations

from typing import Any

from infinite_sets import InfiniteSet, everything
from mdtraj import Trajectory
from narupa.trajectory import FrameData

from narupatools.frame.frame_source import TrajectorySource

from .converter import mdtraj_trajectory_to_frame


class MDTrajTrajectory(TrajectorySource):
    """MDTraj trajectory playback."""

    def __init__(self, trajectory: Trajectory, /):
        super().__init__()
        self._trajectory = trajectory

    def get_frame(  # noqa: D102
        self, *, index: int, fields: InfiniteSet[str]
    ) -> FrameData:
        return mdtraj_trajectory_to_frame(
            self._trajectory, frame_index=index, fields=everything()
        )

    def __len__(self) -> int:
        return len(self._trajectory)

    @classmethod
    def _create_from_object(cls, obj: Any) -> MDTrajTrajectory:
        if isinstance(obj, Trajectory):
            return MDTrajTrajectory(obj)
        raise NotImplementedError
