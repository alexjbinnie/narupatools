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

import time
from typing import Optional

import pytest

from narupatools.core._playable import Playable


class ExamplePlayable(Playable):
    steps: int
    max_steps: Optional[int]

    def __init__(self):
        super().__init__(playback_interval=0.1)
        self.steps = 0
        self.max_steps = None

    def _advance(self) -> bool:
        self.steps += 1
        return self.max_steps is None or self.steps < self.max_steps

    def _restart(self) -> None:
        pass


@pytest.fixture
def playable():
    playable = ExamplePlayable()
    yield playable
    playable.stop(wait=True)


def test_initial_not_playing(playable):
    assert not playable.is_running
    assert not playable.is_playing
    assert not playable.is_paused


def test_step(playable):
    assert playable.steps == 0
    playable.step()
    assert playable.steps == 1
    playable.step()
    assert playable.steps == 2


def test_run(playable):
    playable.max_steps = 10
    playable.run(block=True)
    assert playable.steps == 10


def test_run_is_running(playable):
    playable.max_steps = 100
    playable.run(block=False)
    time.sleep(0.1)
    assert playable.is_running
    assert playable.is_playing
    assert not playable.is_paused


def test_play(playable):
    playable.play()
    time.sleep(1.0)
    assert playable.steps > 0
    playable.stop(wait=True)


def test_pause(playable):
    playable.max_steps = 100
    playable.run(block=False)
    time.sleep(0.1)
    playable.pause()
    assert playable.is_running
    assert not playable.is_playing
    assert playable.is_paused


@pytest.mark.timeout(4)
@pytest.mark.parametrize("interval", [0.1, 0.25])
def test_playback_interval(playable, interval):
    playable.playback_interval = interval
    assert playable.playback_interval == pytest.approx(interval)
    playable.run(block=False)
    time.sleep(2)
    # Accept a 10% error
    assert playable.steps == pytest.approx(2.0 / interval, rel=1e-1)


@pytest.mark.timeout(4)
@pytest.mark.parametrize("rate", [10, 50])
def test_playback_rate(playable, rate):
    playable.playback_rate = rate
    assert playable.playback_rate == pytest.approx(rate)
    playable.run(block=False)
    time.sleep(2)
    # Accept a 20% error
    assert playable.steps == pytest.approx(2.0 * rate, rel=0.5)


def test_play_twice(playable):
    playable.play()
    time.sleep(0.2)
    playable.play()
    assert playable.is_playing


def test_pause_then_play(playable):
    playable.play()
    time.sleep(0.2)
    playable.pause()
    assert playable.is_paused


def test_run_twice(playable):
    playable.run(block=False)
    with pytest.raises(RuntimeError):
        playable.run(block=False)
