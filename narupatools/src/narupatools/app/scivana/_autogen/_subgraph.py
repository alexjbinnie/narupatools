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

import inspect
from typing import Any, ClassVar, Dict, List, Mapping, Type, TypeVar

from narupatools.override import override
from narupatools.state import SharedStateObject
from narupatools.state.typing import Serializable

_TClass = TypeVar("_TClass", bound="SubgraphObject")


class SubgraphObject(SharedStateObject):
    _subgraph_ids: ClassVar[List[str]] = []
    _subclasses: List[Type[SubgraphObject]] = []

    def __init__(self, **kwargs: Serializable):
        super().__init__(**kwargs)

    @classmethod
    @override
    def deserialize(cls: Type[_TClass], value: Serializable) -> _TClass:  # noqa: D102
        if isinstance(value, Mapping):
            subgraph_type = value.get("type", None)
            if subgraph_type is None:
                return super().deserialize(value)
            value = dict(**value)
            del value["type"]
            for subclass in cls._subclasses:
                if subgraph_type in subclass._subgraph_ids:
                    return subclass.deserialize(value)  # type: ignore
            return super().deserialize(value)
        raise ValueError

    @override
    def serialize(self) -> Dict[str, Serializable]:
        serialized = super().serialize()
        serialized["type"] = type(self)._subgraph_ids[0]
        return serialized

    def __init_subclass__(cls, **kwargs: Any) -> None:
        super().__init_subclass__(**kwargs)
        for base in inspect.getmro(cls):
            if (
                base is not cls
                and base is not SubgraphObject
                and issubclass(base, SubgraphObject)
            ):
                base._subclasses = base._subclasses + [cls]
        cls._subclasses = []
