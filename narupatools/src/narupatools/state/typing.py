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

"""Typing for the shared state code."""

from __future__ import annotations

from abc import abstractmethod
from typing import (
    Iterator,
    MutableMapping,
    Protocol,
    Type,
    TypeVar,
    Union,
    runtime_checkable,
)


class SerializableIterable(Protocol):
    """Protocol for an iterable collection of serializable values."""

    def __iter__(self) -> Iterator["Serializable"]:
        ...


class SerializableMapping(Protocol):
    """Protocol for a mapping of string keys to serializable values."""

    def __getitem__(self, key: str) -> "Serializable":
        ...

    def __len__(self) -> int:
        ...

    def __iter__(self) -> Iterator[str]:
        ...


_TClass = TypeVar("_TClass")


@runtime_checkable
class SerializableObject(Protocol):
    """
    Protocol for objects which can be converted to and from a serializable form.

    Types that are serializable by protobuf consist of primitive types, lists and
    dictionaries.
    """

    @classmethod
    @abstractmethod
    def deserialize(cls: Type[_TClass], value: Serializable) -> _TClass:
        """
        Deserialize a given object to an instance of this class.

        Keys which share a name with a property will use that property's setter, whilst
        other keys will be stored as arbitrary data which can be accessed using an item
        accessor.

        :param value: Object to deserialize to the given object.
        """
        raise NotImplementedError

    @abstractmethod
    def serialize(self) -> Serializable:
        """
        Serialize this object to a form that is serializable to protobuf.

        All properties and arbitrary data will be writen to the dictionary.
        """
        raise NotImplementedError


# Define the Serializable type, which is recursively defined as simple primitive types
# and lists/dictionaries of them.
Serializable = Union[
    None,
    str,
    int,
    float,
    bool,
    SerializableIterable,
    SerializableMapping,
    SerializableObject,
]

SerializableDictionary = MutableMapping[str, Serializable]
