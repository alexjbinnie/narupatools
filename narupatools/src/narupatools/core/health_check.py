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

"""Checks background threads to ensure they are still running correctly."""

from abc import ABCMeta
from concurrent.futures import Future
from typing import Optional


class HealthCheck(metaclass=ABCMeta):
    """
    Base class for an object whose health can be checked.

    Health checking is for objects which run tasks in the background, and is used to
    check if there has been an exception on any background threads.
    """

    def health_check(self) -> None:
        """
        Check background tasks to ensure they have not encountered an exception.

        Calling this allows exceptions raised on background threads to be thrown on the
        main thread. It is a useful tool for checking if an object has crashed silently
        in the background.

        It is recommended to call this periodically from the main thread to ensure an
        object is healthy. If no exceptions have occurred, this method will do nothing.
        """
        raise NotImplementedError

    @classmethod
    def check_task(cls, task: Optional[Future]) -> None:
        """
        Check if a task has thrown an exception, and reraise it if so.

        :param task: Background task to check for an exception.
        """
        if task is not None and task.done():
            task.result()
