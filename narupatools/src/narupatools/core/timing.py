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
from typing import Callable


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
