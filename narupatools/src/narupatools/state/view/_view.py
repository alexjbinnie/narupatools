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

"""
A view to a shared state, which does not contain data explicitly.

A shared state dictionary is constantly updating, but it is useful to keep track of
references to individual keys.

A SharedStateView represents a view to a SerializableDictionary, with more useful
utility methods than the simple wrapper SharedStateClientWrapper and
SharedStateServerWrapper. When setting keys, either a normal Serializable value can be
provided, or alternatively any object which implements SerializableObject. This removes
the onus from the user from having to explicitly convert certain objects such as
selections into a serializable form before insertion into the dictionary.

Because the dictionary could be changing, accessing a key directly through the view does
not actually read the value from the dictionary. Instead, it returns a
SharedStateReference, which stores both the key you requested and the dictionary itself.
For more information, see the documentation for SharedStateReference.

A reference to a whole set of keys with a common prefix can be obtained by using the
collection() method.

A snapshot of the current state of the dictionary can be obtained by calling the
snapshot() method. This returns a copy of the dictionary, which is no longer tied to the
original and hence is a snapshot of its value at that point in time.
"""
import abc
from abc import ABC
from contextlib import contextmanager
from typing import (
    AbstractSet,
    Any,
    Dict,
    Generator,
    Generic,
    Iterator,
    Mapping,
    TypeVar,
    Union,
)

from narupatools.state.typing import Serializable, SerializableDictionary

from .._serializable_object import SerializableObject
from ._reference import SharedStateReference

TValue = TypeVar("TValue", bound=Union[SerializableObject, Serializable])


class SharedStateView(ABC, Generic[TValue], Mapping[str, SharedStateReference[TValue]]):
    """
    Represents a view of a shared state dictionary.

    Accessing specific keys will not return their current
    value, but instead a SharedStateReference to a specific key.
    """

    def __init__(self, dictionary: SerializableDictionary, /):
        """
        Create a view of a shared state dictionary, exposing the values as references.

        :param dictionary: Shared state dictionary that implements getting, updating and
                           removing keys.
        """
        self._dictionary = dictionary

    def snapshot(self) -> Dict[str, TValue]:
        """Return a snapshot of each of the current items in this collection."""
        keys = self._keys()
        return {key: self._make_snapshot(key) for key in keys}

    def set(self, key: str, snapshot: TValue, /) -> SharedStateReference[TValue]:
        """Insert a value into the shared state."""
        key = self._resolve_key(key)
        self[key] = snapshot
        return self._make_reference(key)

    def update(
        self, key: str, /, **kwargs: Serializable
    ) -> SharedStateReference[TValue]:
        """
        Update a value by adding the provided key-value pairs.

        This only works if the value supports keyword assignment.

        :param key: Key to modify.
        :param kwargs: Values to add to the value.
        """
        key = self._resolve_key(key)
        self[key].update(**kwargs)
        return self._make_reference(key)

    @contextmanager
    def modify(self, key: str, /) -> Generator[TValue, None, None]:
        """Get the current value, and apply changes after it is modified."""
        key = self._resolve_key(key)
        with self[key].modify() as snapshot:
            yield snapshot

    def __getitem__(self, key: str, /) -> SharedStateReference[TValue]:
        """
        Get a reference to the item with the given key.

        :param key: The key to store in the dictionary
        :raises KeyError: The given key is not present in the shared state dictionary
        """
        key = self._resolve_key(key)
        if not isinstance(key, str):
            raise KeyError

        return self._make_reference(key)

    def __setitem__(self, key: str, value: TValue) -> None:
        """
        Set a shared state value with the given key.

        :param key: The key to alter in the dictionary
        :param value: The value to be inserted into the dictionary. If it implements
                      ``SharedStateSnapshot``, then the ``serialize()`` method is called
                      to generate the actual value to store
        :raises KeyError: The given key is not present in the shared state dictionary
        """
        key = self._resolve_key(key)
        if not isinstance(key, str):
            raise KeyError

        if isinstance(value, SerializableObject):
            self._dictionary[key] = value.serialize()
        else:
            self._dictionary[key] = value  # type: ignore

    def __delitem__(self, key: str) -> None:
        """
        Remove the item with the given key from the shared state dictionary.

        :param key: The key to remove from the dictionary.

        :raises KeyError: The given key is not present in the shared state dictionary.
        """
        key = self._resolve_key(key)
        if not isinstance(key, str):
            raise KeyError

        del self._dictionary[key]

    def __contains__(self, item: Any) -> bool:
        if isinstance(item, str):
            return self._resolve_key(item) in self._dictionary
        elif isinstance(item, SharedStateReference):
            return item.key in self._dictionary
        return False

    def __len__(self) -> int:
        return len(self._keys())

    def _keys(self) -> AbstractSet[str]:
        return self._dictionary.keys()

    def _resolve_key(self, key: str, /) -> str:
        """Resolves a key, so collections can append their prefix if required."""
        return key

    @abc.abstractmethod
    def _make_reference(self, full_key: str, /) -> SharedStateReference[TValue]:
        """
        Create a reference to the item with the given key in this collection.

        :param full_key: The exact key to be referenced. This should begin with the
                         prefix.
        """
        raise NotImplementedError

    def _make_snapshot(self, full_key: str, /) -> TValue:
        """
        Get the current value stored for the given key.

        :param full_key: The exact key to be referenced. This should begin with the
                         prefix.
        """
        return self._make_reference(full_key).snapshot()

    def __repr__(self) -> str:
        return f"<{type(self).__name__} {len(self)} item(s)>"

    def clear(self) -> None:
        """Remove all keys in this dictionary."""
        keys = self._keys()
        for key in list(keys):
            del self[key]

    def __iter__(self) -> Iterator[str]:
        yield from self._keys()
