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

from typing import Iterator, MutableMapping, Protocol, Union


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


# Define the Serializable type, which is recursively defined as simple primitive types
# and lists/dictionaries of them.
Serializable = Union[
    None, str, int, float, bool, SerializableIterable, SerializableMapping
]

SerializableDictionary = MutableMapping[str, Serializable]
