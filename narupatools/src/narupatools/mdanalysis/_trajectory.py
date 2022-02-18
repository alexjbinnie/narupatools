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

"""Adapts a MDAnalysis trajectory for trajectory playback."""

from typing import Any

from infinite_sets import InfiniteSet, everything
from MDAnalysis import Universe
from narupa.trajectory import FrameData

from narupatools.frame import TrajectorySource
from narupatools.mdanalysis import mdanalysis_atomgroup_to_frame


class MDAnalysisTrajectory(TrajectorySource):
    """Wrapper around an MDAnalysis universe to serve as a trajectory source."""

    def __init__(self, universe: Universe):
        self.universe = universe

    def get_frame(  # noqa: D102
        self, *, index: int, fields: InfiniteSet[str] = everything()
    ) -> FrameData:
        _ = self.universe.trajectory[index]
        return mdanalysis_atomgroup_to_frame(self.universe.atoms, fields=fields)

    def __len__(self) -> int:
        return len(self.universe.trajectory)

    @classmethod
    def _create_from_object(cls, obj: Any) -> TrajectorySource:
        if isinstance(obj, Universe):
            return MDAnalysisTrajectory(obj)
        raise NotImplementedError
