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
from narupatools.core.trajectory import TrajectoryPlayback
from narupatools.mdanalysis import mdanalysis_universe_to_frame


class ASETrajectoryPlayback(TrajectoryPlayback):
    """Trajectory playback using one ore more  ASE `Atoms` objects."""

    def __init__(
        self,
        trajectory: Sequence[Atoms],
        *,
        universe: Optional[Universe] = None,
        playback_interval: float = 0.1,
        looping: bool = True,
    ):
        """
        Create a new playback of the given trajectory of ASE `Atoms` objects.

        :param trajectory: Collection of one or more ASE `Atoms` objects representing
                           a trajectory.
        :param universe: Optional MDAnalysis universe with topology information.
        :param playback_interval: The interval between consecutive trajectory frames as
                                  it is played back, in seconds.
        :param looping: Should playback restart from the beginning when the end of the
                        trajectory is reached.
        """
        super().__init__(playback_interval=playback_interval, looping=looping)
        self._trajectory = trajectory
        self._universe = universe

    def _trajectory_length(self) -> int:
        return len(self._trajectory)

    def current_atoms(self) -> Atoms:
        """Currently selected Atoms object in the trajectory."""
        return self._trajectory[self.index]

    def get_frame(self, fields: InfiniteSet[str]) -> FrameData:  # noqa: D102
        frame = FrameData()
        if self._universe:
            ase_atoms_to_frame(self._trajectory[self.index], fields=fields, frame=frame)
            added_fields = set(frame.arrays.keys()) | set(frame.values.keys())
            mdanalysis_universe_to_frame(
                self._universe, fields=fields - added_fields, frame=frame
            )
            return frame
        else:
            ase_atoms_to_frame(self._trajectory[self.index], fields=fields, frame=frame)
        return frame
