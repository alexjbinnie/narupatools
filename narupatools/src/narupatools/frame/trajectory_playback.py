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

"""Base class for trajectories that can be played back."""
from typing import Any, Union

from infinite_sets import InfiniteSet, everything
from narupa.trajectory import FrameData

from narupatools.core.event import Event, EventListener
from narupatools.core.playable import Playable
from narupatools.frame.frame_source import (
    FrameSourceWithNotify,
    OnFieldsChangedCallback,
    TrajectorySource,
)


class TrajectoryPlayback(Playable, FrameSourceWithNotify):
    """Allows trajectories to be played back."""

    _looping: bool
    _on_fields_changed: Event[OnFieldsChangedCallback]

    def __init__(
        self,
        trajectory: Union[TrajectorySource, Any],
        *,
        playback_interval: float = 0.1,
        looping: bool = True,
    ):
        """
        Create a new playback for a trajectory.

        :param playback_interval: The interval between consecutive trajectory frames as
                                  it is played back, in seconds.
        :param looping: Should playback restart from the beginning when the end of the
                        trajectory is reached.
        """
        super().__init__(playback_interval=playback_interval)
        if isinstance(trajectory, TrajectorySource):
            self.trajectory = trajectory
        else:
            traj = TrajectorySource.create_from_object(trajectory)
            if traj is None:
                raise TypeError(f"Failed to convert {trajectory} to a TrajectorySource")
            self.trajectory = traj
        self._looping = looping
        self._index = 0
        self._on_fields_changed = Event(OnFieldsChangedCallback)

    @property
    def on_field_changed(self) -> EventListener[OnFieldsChangedCallback]:  # noqa: D102
        return self._on_fields_changed

    @property
    def looping(self) -> bool:
        """Should playback restart from the beginning when the end is reached?"""
        return self._looping

    @looping.setter
    def looping(self, looping: bool) -> None:
        self._looping = looping

    def get_frame(  # noqa: D102
        self, fields: InfiniteSet[str] = everything()
    ) -> FrameData:
        return self.trajectory.get_frame(index=self.index, fields=fields)

    @property
    def index(self) -> int:
        """Index of the current frame in the trajectory."""
        return self._index

    @index.setter
    def index(self, value: int) -> None:
        if value < 0 or value >= len(self.trajectory):
            raise IndexError(f"Cannot set trajectory index {value}")
        self._index = value
        self._on_fields_changed.invoke(fields=everything())

    def _advance(self) -> bool:
        self._index += 1
        if self._index >= len(self.trajectory):
            if self.looping:
                self._index = 0
                self._on_fields_changed.invoke(fields=everything())
            else:
                self._index = len(self.trajectory) - 1
                self._on_fields_changed.invoke(fields=everything())
                return False
        else:
            self._on_fields_changed.invoke(fields=everything())
        return True

    def _restart(self) -> None:
        self._index = 0
        self._on_fields_changed.invoke(fields=everything())
