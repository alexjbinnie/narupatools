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

"""Protocol for an interactions source."""

from __future__ import annotations

from abc import abstractmethod
from typing import Mapping

from narupa.imd import ParticleInteraction
from typing_extensions import Protocol


class InteractionsSource(Protocol):
    """
    General protocol for an object which can provide interactions.

    Exists to generalize Narupa's ImdStateWrapper so interactions can be provided from
    other sources that are not based on a Narupa server.
    """

    @property
    @abstractmethod
    def active_interactions(self) -> Mapping[str, ParticleInteraction]:
        """Key-indexed set of interactions that should be applied to the system."""
        ...


def wrap_interaction_source(
    dict: Mapping[str, ParticleInteraction], /
) -> InteractionsSource:
    """
    Wrap a pythonic dictionary as an interaction.

    :param dict: Generic mapping of string keys to interactions.
    """
    return _InteractionsSourceWrapper(dict)


class _InteractionsSourceWrapper(InteractionsSource):
    def __init__(self, wrapped: Mapping[str, ParticleInteraction]):
        self._wrapped = wrapped

    @property
    def active_interactions(self) -> Mapping[str, ParticleInteraction]:  # noqa:D102
        return self._wrapped
