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

"""Wrapper around an MDAnalysis universe."""

from infinite_sets import InfiniteSet
from MDAnalysis import Universe
from narupa.trajectory import FrameData

from narupatools.frame._frame_source import FrameSource
from narupatools.mdanalysis import mdanalysis_universe_to_frame


class MDAnalysisSystem(FrameSource):
    """MDAnalysis universe that can be broadcast on a session."""

    def __init__(self, universe: Universe):
        self._universe = universe

    def get_frame(self, fields: InfiniteSet[str]) -> FrameData:  # noqa: D102
        return mdanalysis_universe_to_frame(self._universe)
