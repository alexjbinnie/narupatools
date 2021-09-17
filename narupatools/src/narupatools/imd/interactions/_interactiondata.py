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

from __future__ import annotations

from typing import Dict, Mapping, Type

from narupatools.core.properties import bool_property, numpy_property, str_property, float_property
from narupatools.override import override
from narupatools.state import SharedStateObject
from narupatools.state.typing import Serializable

_InteractionData_Types: Dict[str, Type[InteractionData]] = {}
_InteractionFeedback_Types: Dict[str, Type[InteractionFeedback]] = {}


class InteractionData(SharedStateObject):
    """
    Interaction data that can be serialized to a shared state.

    This class is dynamically subclassed based on the interaction_type field.
    """

    @numpy_property(dtype=int)
    def particles(self) -> None:
        """Indices of the particles affected by this interaction."""
        ...

    @str_property
    def interaction_type(self) -> str:
        """Internal type of the interaction."""
        ...

    @bool_property
    def reset_velocities(self) -> bool:
        """Should this interaction perform velocity reset afterwards?"""
        ...

    @bool_property
    def mass_weighted(self) -> bool:
        """Is this interaction mass-weighted?"""
        ...

    @staticmethod
    def register_interaction_type(
            interaction_type: str, python_type: Type[InteractionData], /
    ) -> None:
        """Register a new subclass of InteractionData and the interaction_type it affects."""
        _InteractionData_Types[interaction_type] = python_type

    @classmethod
    @override
    def deserialize(cls, value: Serializable) -> InteractionData:  # noqa: D102
        if cls is InteractionData and isinstance(value, Mapping):
            interaction_type = value["interaction_type"]
            return _InteractionData_Types[interaction_type].deserialize(value)
        return super().deserialize(value)  # type: ignore[return-value]


class InteractionFeedback(SharedStateObject):
    """
    Data which is produced about the interaction and streamed back to clients.

    This class is dynamically subclassed based on the interaction_type field.
    """


    @str_property
    def interaction_type(self) -> str:
        """Internal type of the interaction."""
        ...


    @float_property
    def work(self) -> str:
        """Internal type of the interaction."""
        ...

    @float_property
    def potential_energy(self) -> bool:
        """Should this interaction perform velocity reset afterwards?"""
        ...

    @staticmethod
    def register_interaction_type(
            interaction_type: str, python_type: Type[InteractionFeedback], /
    ) -> None:
        """Register a new subclass of InteractionData and the interaction_type it affects."""
        _InteractionFeedback_Types[interaction_type] = python_type

    @classmethod
    @override
    def deserialize(cls, value: Serializable) -> InteractionData:  # noqa: D102
        if cls is InteractionFeedback and isinstance(value, Mapping):
            interaction_type = value["interaction_type"]
            return _InteractionFeedback_Types[interaction_type].deserialize(value)
        return super().deserialize(value)  # type: ignore[return-value]
