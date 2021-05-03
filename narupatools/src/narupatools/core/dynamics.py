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

"""Base implementation of simulation dynamics."""

from __future__ import annotations

from abc import ABCMeta, abstractmethod
from concurrent.futures import Future
from typing import Optional, Protocol, Union

from infinite_sets import InfiniteSet
from narupa.trajectory import FrameData

from narupatools.core.event import Event, EventListener
from .playable import Playable
from ..frame import (
    SimulationElapsedSteps,
    SimulationElapsedTime,
    SimulationTotalSteps,
    SimulationTotalTime,
)
from ..frame.frame_source import FrameSource


class OnResetCallback(Protocol):
    """Callback for when a simulation is reset."""

    def __call__(self) -> None:  # noqa: D102
        ...


class OnPreStepCallback(Protocol):
    """Callback for just before an MD step is taken."""

    def __call__(self) -> None:  # noqa: D102
        ...


class OnPostStepCallback(Protocol):
    """Callback for just after an MD step is taken."""

    def __call__(self) -> None:  # noqa: D102
        ...


class SimulationDynamics(Playable, FrameSource, metaclass=ABCMeta):
    """
    Base class for an implementation of dynamics driven by Narupa.

    This implements all the common commands for controlling a simulation, as well as
    running the simulation in both blocking and non-blocking modes. Subclasses of this
    only have to override _step_internal and _reset_internal to implement stepping
    forward a single step in the simulation and reseting the entire simulation back to
    its initial state.
    """

    _remaining_steps: Optional[int]
    _previous_total_steps: int
    _previous_total_time: float
    _on_reset: Event[OnResetCallback]
    _on_pre_step: Event[OnPreStepCallback]
    _on_post_step: Event[OnPostStepCallback]

    def __init__(self, *, playback_interval: float):
        """
        Create a SimulationDynamics object.

        Should be called by subclasses, and not called directly.

        :param playback_interval: Interval at which dynamics will be run, in seconds.
        """
        super().__init__(playback_interval=playback_interval)
        self._on_reset = Event(OnResetCallback)
        self._on_pre_step = Event(OnPreStepCallback)
        self._on_post_step = Event(OnPostStepCallback)
        self._previous_total_time = 0.0
        self._previous_total_steps = 0
        self._elapsed_time = 0.0
        self._elapsed_steps = 0
        self._remaining_steps = None

    @property
    def on_reset(self) -> EventListener[OnResetCallback]:
        """
        Event triggered when dynamics is reset.

        This can be used when the behaviour of reset that the current dynamics performs
        is insufficient, or the user would like to perform a more detailed reset such as
        reinitializing velocities.
        """
        return self._on_reset

    @property
    def on_pre_step(self) -> EventListener[OnPreStepCallback]:
        """Event triggered before each step of the dynamics is run."""
        return self._on_pre_step

    @property
    def on_post_step(self) -> EventListener[OnPostStepCallback]:
        """Event triggered after each step of the dynamics is run."""
        return self._on_post_step

    def reset(self) -> None:
        """
        Reset the simulation.

        The behaviour of this depends on the implementation of the dynamics,
        but commonly an implementation of `SimulationDynamics` should record the initial
        state of the simulation on initialization, and use this to reset the simulation.

        This method is called whenever a client runs the reset command,
        described in :ref:`narupa.trajectory.frame_server`.
        """
        super().restart()

    def _restart(self) -> None:
        self._previous_total_time += self.elapsed_time
        self._previous_total_steps += self.elapsed_steps
        self._elapsed_time = 0
        self._elapsed_steps = 0
        self._reset_internal()
        self._on_reset.invoke()

    def run(  # type: ignore
        self, steps: Optional[int] = None, block: Optional[bool] = None
    ) -> Union[bool, Future[bool]]:
        """
        Run the dynamics.

        :raises ValueError: The number of steps was negative or 0.
        :param steps: The number of steps to run for, or None to run indefinitely. If
                      provided, must be larger than 0.
        :param block: Should this be run in this thread (block=True) or in a background
                      thread (block=False)
        :return: If run in blocking mode, returns True if the dynamics completed and
                 False if it was interrupted. If not run in blocking mode, returns a
                 Future with the same result.
        """
        if steps is not None and steps <= 0:
            raise ValueError("Cannot run for less than 1 steps.")
        if block is None:
            block = steps is not None
        self._remaining_steps = steps
        return super().run(block)

    def _advance(self) -> bool:
        self._on_pre_step.invoke()
        self._step_internal()
        self._elapsed_steps += 1
        self._elapsed_time += self.time_step
        self._on_post_step.invoke()
        if self._remaining_steps is not None:
            self._remaining_steps -= 1
            return self._remaining_steps > 0
        else:
            return True

    @abstractmethod
    def _step_internal(self) -> None:
        """Step the dynamics forward by a single timestep."""
        raise NotImplementedError()

    @abstractmethod
    def _reset_internal(self) -> None:
        """Reset the simulation to its initial state."""
        raise NotImplementedError()

    @property
    def elapsed_time(self) -> float:
        """
        Elapsed time of the simulation since initialization/last reset in picoseconds.

        When the system is reset, this value is reset to 0.
        """
        return self._elapsed_time

    @property
    def elapsed_steps(self) -> int:
        """
        Elapsed number of steps of the simulation since initialization/last reset.

        When the system is reset, this value is reset to 0.
        """
        return self._elapsed_steps

    @property
    def total_time(self) -> float:
        """
        Total time of the simulation in picoseconds.

        This includes all times the simulation has been reset.
        """
        return self._previous_total_time + self.elapsed_time

    @property
    def total_steps(self) -> int:
        """
        Total number of steps of the simulation.

        This includes all times the simulation has been reset.
        """
        return self._previous_total_steps + self.elapsed_steps

    @property
    @abstractmethod
    def time_step(self) -> float:
        """
        Current time step of the simulation in picoseconds.

        :raises AttributeError: Cannot get the current time step for this dynamics.
        """
        raise AttributeError()

    @property
    def temperature(self) -> float:
        """
        Current temperature of the dynamics in Kelvin.

        :raises AttributeError: Temperature is not defined for this dynamics.
        """
        raise AttributeError()

    @abstractmethod
    def _get_frame(self, fields: InfiniteSet[str]) -> FrameData:
        pass

    def get_frame(self, fields: InfiniteSet[str]) -> FrameData:
        """
        Get the current state of the system as a Narupa `FrameData`.

        :param fields: Collection of keys to include in the FrameData
        :return: Narupa `FrameData` populated with requested fields.
        """
        frame = self._get_frame(fields)
        SimulationElapsedTime.set(frame, self.elapsed_time)
        SimulationElapsedSteps.set(frame, self.elapsed_steps)
        SimulationTotalTime.set(frame, self.total_time)
        SimulationTotalSteps.set(frame, self.total_steps)
        return frame
