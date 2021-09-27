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

"""Code for handling FrameData and other related objects."""

from ._converter import FrameConverter, convert, frame_to_pdb_string
from ._frame import NarupaFrame
from ._frame_producer import FrameProducer
from ._frame_source import (
    FrameSource,
    FrameSourceWithNotify,
    OnFieldsChangedCallback,
    TrajectorySource,
)
from ._patch import *  # noqa: F401, F403
from ._trajectory_playback import TrajectoryPlayback

__all__ = [
    "convert",
    "NarupaFrame",
    "TrajectoryPlayback",
    "FrameSource",
    "FrameSourceWithNotify",
    "OnFieldsChangedCallback",
    "TrajectorySource",
    "FrameConverter",
    "FrameProducer",
    "frame_to_pdb_string",
]
