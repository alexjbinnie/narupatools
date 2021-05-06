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

import pytest
from infinite_sets import InfiniteSet
from narupa.trajectory import FrameData

from narupatools.core.dynamics import SimulationDynamics
from testing import assert_event_called, assert_event_not_called


class ExampleDynamics(SimulationDynamics):
    def _get_frame(self, fields: InfiniteSet[str]) -> FrameData:
        return FrameData()

    def __init__(self, *, playback_interval: float = 0.0):
        super().__init__(playback_interval=playback_interval)

    def _reset_internal(self):
        pass

    def _step_internal(self):
        pass

    @property
    def timestep(self) -> float:
        return 1.0


@pytest.fixture
def dynamics():
    dynamics = ExampleDynamics()
    yield dynamics
    dynamics.stop(wait=True)


def test_initial_not_running(dynamics):
    assert not dynamics.is_running


def test_initial_elapsed_time(dynamics):
    assert dynamics.elapsed_time == 0


def test_initial_total_time(dynamics):
    assert dynamics.total_time == 0


def test_initial_elapsed_steps(dynamics):
    assert dynamics.elapsed_steps == 0


def test_initial_total_steps(dynamics):
    assert dynamics.total_steps == 0


@pytest.mark.parametrize("steps", [1, 5, 10, 25])
def test_run_elapsed_and_total_steps(dynamics, steps):
    dynamics.run(steps)
    assert dynamics.elapsed_steps == steps
    assert dynamics.total_steps == steps


@pytest.mark.parametrize("steps", [1, 5, 10, 25])
def test_run_elapsed_and_total_time(dynamics, steps):
    dynamics.run(steps)
    assert dynamics.elapsed_time == pytest.approx(steps)
    assert dynamics.total_time == pytest.approx(steps)


@pytest.mark.parametrize("steps", [1, 5, 10, 25])
def test_reset_elapsed_and_total_steps(dynamics, steps):
    dynamics.run(steps)
    dynamics.reset()
    assert dynamics.elapsed_steps == 0
    assert dynamics.total_steps == steps
    dynamics.run(steps)
    assert dynamics.elapsed_steps == steps
    assert dynamics.total_steps == 2 * steps


@pytest.mark.parametrize("steps", [1, 5, 10, 25])
def test_reset_elapsed_and_total_time(dynamics, steps):
    dynamics.run(steps)
    dynamics.reset()
    assert dynamics.elapsed_steps == pytest.approx(0.0)
    assert dynamics.total_steps == pytest.approx(steps)
    dynamics.run(steps)
    assert dynamics.elapsed_steps == pytest.approx(steps)
    assert dynamics.total_steps == pytest.approx(2 * steps)


def test_cancel_run_non_blocking(dynamics):
    dynamics.run(block=False)
    assert dynamics.is_running
    dynamics.stop(wait=True)
    assert dynamics.is_running is False


def test_cancel_is_atomic(dynamics):
    dynamics.stop()
    dynamics.stop()


def test_run_twice(dynamics):
    dynamics.run(block=False)
    with pytest.raises(RuntimeError):
        dynamics.run(block=False)
    dynamics.stop(wait=True)


def test_play(dynamics):
    dynamics.run(block=False)
    assert dynamics.is_playing
    dynamics.pause()
    assert dynamics.is_playing is False
    dynamics.play()
    assert dynamics.is_playing
    dynamics.stop(wait=True)


def test_play_twice(dynamics):
    dynamics.play()
    assert dynamics.is_playing
    dynamics.play()
    assert dynamics.is_playing


def test_step(dynamics):
    assert dynamics.total_steps == 0
    dynamics.step()
    assert dynamics.total_steps == 1


@pytest.mark.parametrize("nsteps", [1, 5, 10, 25])
def test_multiple_steps(dynamics, nsteps):
    assert dynamics.total_steps == 0
    for _ in range(nsteps):
        dynamics.step()
    assert dynamics.total_steps == nsteps


def test_step_pauses_simulation(dynamics):
    dynamics.run(block=False)
    time.sleep(0.1)
    assert dynamics.total_steps > 0
    dynamics.step()
    assert dynamics.is_paused


def test_reset(dynamics):
    with assert_event_not_called(dynamics.on_reset):
        dynamics.run(10)

    with assert_event_called(dynamics.on_reset):
        dynamics.reset()
