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
View on a subset of keys in a shared state.

A SharedStateCollectionView is a view on a subset of keys belonging in a dictionary. See
the documentation for SharedStateDictionaryView for more information.

A collection view is defined by a prefix such as 'selection.', which indicates what
subset of keys it encapsulates. All methods referring to keys are prepended by this
prefix if not so already, so accessing collection['abc'] for a collection with the
prefix 'selection.' returns a reference to the key 'selection.abc'

A collection, like a reference, can also specify that snapshots of its items should be
converted to a specific class that derives from SerializableObject. See
SharedStateReference for more information.
"""
from __future__ import annotations

import uuid
from typing import AbstractSet, Generic, Optional, Type, TypeVar, Union

from narupatools.state.typing import Serializable, SerializableDictionary

from ..serializable_object import SerializableObject
from .reference import SharedStateReference
from .view import SharedStateView

TValue = TypeVar("TValue", bound=Union[SerializableObject, Serializable])

TSpecificType = TypeVar("TSpecificType", bound=SerializableObject)


class SharedStateCollectionView(SharedStateView[TValue], Generic[TValue]):
    """View on a collection of items in a shared state with a given prefix."""

    def __init__(
        self,
        view: SerializableDictionary,
        prefix: str,
        snapshot_type: Optional[Type[SerializableObject]] = None,
    ):
        """
        Create a reference to a collection of items with a given prefix.

        :param view: A view of the shared state dictionary that handles getting,
                     updating and removing keys
        :param prefix: The prefix that keys must have to be considered in this
                       collection. For example, a prefix of 'selection.' matches all
                       keys of the form 'selection.*'
        :param snapshot_type: The class that snapshots of its in this collection should
                              implement
        """
        super().__init__(view)
        self.prefix = prefix
        self._snapshot_class = snapshot_type

    @classmethod
    def untyped_collection(
        cls, view: SerializableDictionary, prefix: str
    ) -> "SharedStateCollectionView[Serializable]":
        """Create a view of a set of keys that does not assume any specific type."""
        return SharedStateCollectionView(view, prefix)

    @classmethod
    def typed_collection(
        cls,
        view: SerializableDictionary,
        prefix: str,
        snapshot_class: Type[TSpecificType],
    ) -> "SharedStateCollectionView[TSpecificType]":
        """Create a view of a set of keys with values of the specified type."""
        return SharedStateCollectionView(view, prefix, snapshot_class)

    def set(
        self, key: str, snapshot: Optional[TValue] = None, /, **kwargs: Serializable
    ) -> SharedStateReference[TValue]:
        """
        Insert a value into the shared state dictionary with the given key.

        This returns a reference to the item that was inserted into the shared state.

        :param key: Key to modify.
        :param snapshot: Value to set.
        :param kwargs: Arbitrary key-value pairs to set on the value.
        :raises ValueError: Can only use kwargs if this is a typed view.
        """
        key = self._resolve_key(key)
        if snapshot is None:
            if self._snapshot_class is None:
                raise ValueError(
                    "Cannot use keywords to instantiate shared state object if "
                    "snapshot_class is not specified"
                )
            obj = self._snapshot_class()
            for keyword, value in kwargs.items():
                obj[keyword] = value  # type: ignore
            return super().set(key, obj)  # type: ignore
        else:
            return super().set(key, snapshot)

    def add(
        self, snapshot: Optional[TValue] = None, /, **kwargs: Serializable
    ) -> SharedStateReference[TValue]:
        """Insert a value into the shared state with an autogenerated key."""
        key = self._resolve_key(self._create_key())
        return self.set(key, snapshot, **kwargs)

    def _keys(self) -> AbstractSet[str]:
        return {key for key in self._dictionary if key.startswith(self.prefix)}

    def _resolve_key(self, key: str) -> str:
        if isinstance(key, str) and not key.startswith(self.prefix):
            return self.prefix + key
        return key

    def _create_key(self) -> str:
        """Generate a unique key in the collection."""
        return self.prefix + str(uuid.uuid4())

    def _make_reference(self, full_key: str) -> SharedStateReference[TValue]:
        """
        Create a reference to a specific item with the given key.

        :param full_key: The exact key to be referenced. This should begin with the
                         prefix.
        """
        if self._snapshot_class is not None:
            return SharedStateReference.typed_reference(
                self._dictionary, full_key, self._snapshot_class  # type: ignore
            )
        else:
            return SharedStateReference.untyped_reference(
                self._dictionary, full_key
            )  # type: ignore
