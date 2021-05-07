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

"""Code for managing timings."""
import time
from contextlib import contextmanager
from typing import Callable

from narupatools.core.event import EventListener


def wait_for(
    predicate: Callable[..., bool], *, timeout: float = 2, check_interval: float = 0.1
) -> None:
    """
    Wait for a given predicate to evaluate to true.

    :param predicate: Event to wait for.
    :param timeout: Time to wait in seconds.
    :param check_interval: Interval to check predicate in seconds.
    :raises TimeoutError: Timeout is reached and the predicate hasn't returned true.
    """
    end_time = time.monotonic() + timeout
    while not (passed := predicate()) and time.monotonic() < end_time:
        time.sleep(check_interval)
    if not passed:
        raise TimeoutError(f"{predicate} did not return True within timeout.")


def wait_for_event(
    event: EventListener, *, timeout: float = 2, check_interval: float = 0.1
):
    """
    Wait for a given event to be invoked.

    :param event: Event to wait for.
    :param timeout: Time to wait in seconds.
    :param check_interval: Interval to check event in seconds.
    :raises TimeoutError: Timeout is reached and the event hasn't been called.
    """
    was_called = False

    def callback(**kwargs):
        nonlocal was_called
        print("CALLBACK!")
        was_called = True

    event.add_callback(callback)
    try:
        wait_for(lambda: was_called, timeout=timeout, check_interval=check_interval)
    except TimeoutError:
        raise TimeoutError(f"{event} was not called within timeout")
    finally:
        event.remove_callback(callback)


@contextmanager
def wait_until_event(
    event: EventListener, *, timeout: float = 2, check_interval: float = 0.1
):
    """
    Wait for a given event to be invoked, as a context manager.

    This is for when you want to run some code followed by waiting for an event, but
    using wait_for_event might miss the event (as the callback is added after you run
    the code). This adds the callback before you run the code.

    :param event: Event to wait for.
    :param timeout: Time to wait in seconds.
    :param check_interval: Interval to check event in seconds.
    :raises TimeoutError: Timeout is reached and the event hasn't been called.
    """
    was_called = False

    def callback(**kwargs):
        nonlocal was_called
        was_called = True

    event.add_callback(callback)

    try:
        was_called = False
        yield
        wait_for(lambda: was_called, timeout=timeout, check_interval=check_interval)
    except TimeoutError:
        print(was_called)
        raise TimeoutError(f"{event} was not called within timeout")
    finally:
        event.remove_callback(callback)
