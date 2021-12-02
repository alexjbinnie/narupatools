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

import pytest
from infinite_sets import InfiniteSet, everything
from narupa.trajectory import FrameData

from narupatools.frame import TrajectoryPlayback, TrajectorySource

ALL_FIELDS = everything()


class TestTrajectory(TrajectorySource):
    def __init__(self):
        super().__init__()

    def __len__(self) -> int:
        return 5

    def get_frame(self, *, index: int, fields: InfiniteSet[str]) -> FrameData:
        data = FrameData()
        data.values["index"] = index
        return data


@pytest.fixture
def trajectory():
    return TestTrajectory()


@pytest.fixture
def playback(trajectory):
    return TrajectoryPlayback(trajectory)


def test_initial_index(playback):
    assert playback.index == 0
    assert playback.get_frame().values["index"] == 0


def test_set_index(playback):
    playback.index = 3
    assert playback.index == 3
    assert playback.get_frame().values["index"] == 3


def test_set_index_out_of_range(playback):
    with pytest.raises(IndexError):
        playback.index = 14


def test_set_index_negative(playback):
    with pytest.raises(IndexError):
        playback.index = -11


def test_increment_step(playback):
    playback.step()
    assert playback.index == 1
    assert playback.get_frame().values["index"] == 1


def test_increment_step_multiple(playback):
    for _ in range(4):
        playback.step()
    assert playback.index == 4
    assert playback.get_frame().values["index"] == 4


def test_looping(playback):
    for _ in range(7):
        playback.step()
    assert playback.index == 2
    assert playback.get_frame().values["index"] == 2


def test_non_looping(playback):
    playback.looping = False
    for _ in range(7):
        playback.step()
    assert playback.index == 4
    assert playback.get_frame().values["index"] == 4


def test_restart(playback):
    for _ in range(4):
        playback.step()
    playback.restart()
    assert playback.index == 0
    assert playback.get_frame().values["index"] == 0
