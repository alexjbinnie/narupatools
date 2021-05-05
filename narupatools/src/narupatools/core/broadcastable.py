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

"""Base class for objects with callbacks when added or removed from a session."""

from __future__ import annotations

from abc import ABCMeta, abstractmethod
from typing import Protocol


class Broadcaster(Protocol):
    """Generalization of an object that can share something such as a simulation."""

    pass


class Broadcastable(metaclass=ABCMeta):
    """Base class for objects that can be broadcast in a session."""

    @abstractmethod
    def start_broadcast(self, broadcaster: Broadcaster) -> None:
        """
        Called when this object is about to be broadcast.

        :param broadcaster: Broadcaster that is about to broadcast this object.
        """
        ...

    @abstractmethod
    def end_broadcast(self, broadcaster: Broadcaster) -> None:
        """
        Called when this object is stopped being broadcasted.

        :param broadcaster: Broadcaster that is about to stop broadcasting this object.
        """
        ...
