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

from narupatools.core.properties import bool_property, numpy_property, str_property
from narupatools.override import override
from narupatools.state import SharedStateObject
from narupatools.state.typing import Serializable

_InteractionParameters_Types: Dict[str, Type[InteractionParameters]] = {}


class InteractionParameters(SharedStateObject):
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
        interaction_type: str, python_type: Type[InteractionParameters], /
    ) -> None:
        """Register a new subclass of InteractionData and the interaction_type it affects."""
        _InteractionParameters_Types[interaction_type] = python_type

    @classmethod
    @override
    def deserialize(cls, value: Serializable) -> InteractionParameters:  # noqa: D102
        if cls is InteractionParameters and isinstance(value, Mapping):
            interaction_type = value["interaction_type"]
            return _InteractionParameters_Types[interaction_type].deserialize(value)
        return super().deserialize(value)  # type: ignore[return-value]
