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

from abc import ABCMeta, abstractmethod
from typing import Dict, Generic, Set, TypeVar

from infinite_sets import InfiniteSet, everything

from narupatools.core.dynamics import SimulationDynamics
from narupatools.physics.typing import Vector3

from ._feature import InteractionFeature

_TDynamics = TypeVar("_TDynamics", bound=SimulationDynamics)


class SetAndClearInteractionFeature(
    Generic[_TDynamics],
    InteractionFeature[_TDynamics],
    metaclass=ABCMeta,
):
    """IMD manager for dynamics that support arbitrary forces on each atom."""

    def __init__(self, dynamics: _TDynamics):
        super().__init__(dynamics)
        self.current_indices: Set[int] = set()

    def _calculate_and_apply_interactions(self) -> None:
        indices: Set[int] = set()
        forces: Dict[int, Vector3] = {}
        for interaction in self.current_interactions.values():
            for index, force in zip(interaction.particle_indices, interaction.forces):
                indices.add(index)
                if index in forces:
                    forces[index] += force
                else:
                    forces[index] = force
        to_clear = self.current_indices - indices
        if to_clear:
            self._clear_forces(to_clear)
        if forces:
            self._set_forces(forces)
        self.current_indices = indices

    @abstractmethod
    def _set_forces(self, forces: Dict[int, Vector3], /) -> None:
        pass

    @abstractmethod
    def _clear_forces(self, indices: InfiniteSet[int] = everything(), /) -> None:
        pass
