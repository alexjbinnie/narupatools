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

"""Base implementation for something that can be played back."""

from __future__ import annotations

import time
from abc import ABCMeta, abstractmethod
from concurrent.futures import Future, ThreadPoolExecutor
from threading import Lock
from typing import Generator, Literal, Optional, Union, overload

from .health_check import HealthCheck


class PlayableAlreadyRunningError(RuntimeError):
    """Error raised when a playable is attempted to be played more than once."""

    def __init__(self) -> None:
        super().__init__(
            "Already running on another thread. Call `stop()` first to "
            "ensure this is not already running"
        )


class Playable(HealthCheck, metaclass=ABCMeta):
    """
    Object that can be run at a target frame rate with playback controls..

    Playables can therefore both represent trajectories (where there are a set number of
    frames which are played through one by one) and simulations (which generate frames
    one at a time).

    Playables support pausing, restarting and stepping, as well as execution in both
    background and foreground threads.
    """

    _playback_interval: float
    _advance_lock: Lock
    _is_running_lock: Lock
    _threads: ThreadPoolExecutor

    _cancelled: bool
    _paused: bool

    _run_task: Optional[Future[bool]]

    def __init__(self, *, playback_interval: float = 0.1):
        """
        Create a new playable.

        :param playback_interval: Target interval to run this playable at.
        """
        self._playback_interval = playback_interval
        self._advance_lock = Lock()
        self._is_running_lock = Lock()
        self._threads = ThreadPoolExecutor(
            max_workers=1, thread_name_prefix=f"{type(self).__name__}_Run"
        )
        self._cancelled = False
        self._paused = False
        self._run_task = None

    @property
    def playback_interval(self) -> float:
        """Time in seconds between individual steps of this playable."""
        return self._playback_interval

    @playback_interval.setter
    def playback_interval(self, value: float) -> None:
        self._playback_interval = value

    @property
    def playback_rate(self) -> float:
        """Number of steps this playable goes through in 1 second."""
        return 1.0 / self._playback_interval

    @playback_rate.setter
    def playback_rate(self, value: float) -> None:
        self._playback_interval = 1.0 / value

    def _yield_at_framerate(self) -> Generator[float, None, None]:
        """Yields at intervals specified by the playback interval."""
        last_yield = time.monotonic() - self.playback_interval
        while True:
            time_since_yield = time.monotonic() - last_yield
            wait_duration = max(0.0, self.playback_interval - time_since_yield)
            time.sleep(wait_duration)
            yield time.monotonic() - last_yield
            last_yield = time.monotonic()

    @property
    def is_running(self) -> bool:
        """
        Is this currently running?

        This is True even if the playback is paused.
        """
        return self._is_running_lock.locked()

    @property
    def is_paused(self) -> bool:
        """Is this currently running but paused?"""
        return self.is_running and self._paused

    @property
    def is_playing(self) -> bool:
        """Is this currently playing and not paused?"""
        return self.is_running and not self._paused

    def restart(self) -> None:
        """Restart this back to its initial state."""
        with self._advance_lock:
            self._restart()

    def play(self) -> None:
        """
        Play the given object.

        If this is not running, run this on a background thread. Else, unpause if this
        is paused.
        """
        if self.is_running:
            self._paused = False
        else:
            self.run(block=False)

    def pause(self, wait: bool = True) -> None:
        """
        Pause this object if it is running.

        This will not cancel the current run if present, merely suspend it until a
        subsequent call to play() is made.

        :param wait: Wait until the pause has actually occurred.
        """
        if self.is_running:
            if wait:
                with self._advance_lock:
                    self._paused = True
            else:
                self._paused = True

    def stop(self, wait: bool = True) -> None:
        """
        Stop the object if it is running.

        If this is running, this object will stop executing at the next available point.

        :param wait: Should this call wait until this object has stopped running before
                     returning?
        """
        if not self.is_running:
            return

        self._cancelled = True

        if wait:
            with self._is_running_lock:
                pass

    def step(self) -> None:
        """Step forward by one step, and then pause the simulation."""
        with self._advance_lock:
            self._advance()
            self._paused = True

    def start(self, block: bool = True) -> None:
        """
        An alias for run().

        :param block: Should this block?
        """
        self.run(block)

    @overload
    def run(self, block: Literal[True]) -> bool:
        ...

    @overload
    def run(self, block: Literal[False]) -> Future[bool]:
        ...

    @overload
    def run(self, block: bool) -> Union[bool, Future[bool]]:
        ...

    def run(self, block: bool = True) -> Union[bool, Future[bool]]:
        """
        Run this object, potentially in a background thread.

        If block is True (as is the default), this object will run in the calling
        thread. If block is False, it will be run in a background thread.

        :param block: Should this block?
        :raises PlayableAlreadyRunningError: Playable is running on another thread.
        :returns: If run in blocking mode, returns True if the playable completed and
                  False if it was stopped. If run in non-blocking mode, returns a Future
                  with the same result.
        """
        if self.is_running:
            raise PlayableAlreadyRunningError()
        if block:
            return self._run()
        else:
            self._run_task = self._threads.submit(self._run)
            return self._run_task

    def _run(self) -> bool:
        with self._is_running_lock:
            self._start_run()
            was_cancelled = False
            for _ in self._yield_at_framerate():
                with self._advance_lock:
                    if self._cancelled:
                        was_cancelled = True
                        break
                    if self._paused:
                        continue
                    if not self._advance():
                        break
            self._end_run()
            self._paused = False
            self._cancelled = False
            return not was_cancelled

    def _start_run(self) -> None:
        pass

    def _end_run(self) -> None:
        pass

    @abstractmethod
    def _advance(self) -> bool:
        ...

    @abstractmethod
    def _restart(self) -> None:
        ...

    def health_check(self) -> None:  # noqa: D102
        HealthCheck.check_task(self._run_task)
