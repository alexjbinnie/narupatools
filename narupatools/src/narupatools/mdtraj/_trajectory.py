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

from mdtraj import Trajectory

from narupatools.core.trajectory import TrajectoryPlayback


class MDTrajTrajectoryPlayback(TrajectoryPlayback):
    """MDTraj trajectory playback."""

    def __init__(self, trajectory: Trajectory, /):
        super().__init__()
        self._trajectory = trajectory

    def _trajectory_length(self) -> int:
        return len(self._trajectory)
