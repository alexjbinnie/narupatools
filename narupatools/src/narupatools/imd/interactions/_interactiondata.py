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

from narupatools.core.properties import numpy_property, str_property, bool_property
from narupatools.state import SharedStateObject
from narupatools.state.typing import Serializable


class InteractionData(SharedStateObject):
    """
    Interaction data that can be serialized to a shared state.

    This class is dynamically subclassed based on the interaction_type field.
    """

    _types: Dict[str, Type[InteractionData]] = {}

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

    @classmethod
    def register_interaction_type(
        cls, interaction_type: str, type: Type[InteractionData]
    ) -> None:
        """Register a new subclass of InteractionData and the interaction_type it affects."""
        cls._types[interaction_type] = type

    @classmethod
    def _get_type(cls, interaction_type: str) -> Type[InteractionData]:
        return cls._types[interaction_type]

    @classmethod
    def deserialize(cls, value: Serializable) -> InteractionData:  # noqa: D102
        if cls is InteractionData and isinstance(value, Mapping):
            interaction_type = value["interaction_type"]
            return cls._get_type(interaction_type).deserialize(value)
        return super().deserialize(value)  # type: ignore[return-value]
