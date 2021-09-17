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
View of a full shared state dictionary.

A SharedStateDictionaryView is a kind of SharedStateView which can generate collections
using the collection() method This specification of SharedStateView into
SharedStateDictionaryView and SharedStateCollectionView ensures that collections
themselves cannot spawn more collection, as well as avoiding the class referencing
itself.
"""

from typing import Optional, Type, TypeVar, overload

from narupatools.state.typing import Serializable

from .._serializable_object import SerializableObject
from ._collection_view import SharedStateCollectionView
from ._reference import SharedStateReference
from ._view import SharedStateView
from ...override import override

TSerializableObjectType = TypeVar("TSerializableObjectType", bound=SerializableObject)


class SharedStateDictionaryView(SharedStateView[Serializable]):
    """Represents a view of a shared state dictionary."""

    @overload
    def collection(
        self, prefix: str, /
    ) -> SharedStateCollectionView[Serializable]:  # noqa: D102
        pass

    @overload
    def collection(
        self,
        prefix: str,
        snapshot_type: Type[TSerializableObjectType],
        /
        # noqa: D102
    ) -> SharedStateCollectionView[TSerializableObjectType]:
        pass

    def collection(
        self, prefix: str, snapshot_type: Optional[Type[SerializableObject]] = None, /
    ) -> SharedStateCollectionView:
        """Get a collection view for all keys with the given prefix."""
        if snapshot_type is not None:
            return SharedStateCollectionView.typed_collection(
                self._dictionary, prefix, snapshot_type
            )
        else:
            return SharedStateCollectionView.untyped_collection(
                self._dictionary, prefix
            )

    @override
    def _make_reference(self, full_key: str, /) -> SharedStateReference[Serializable]:
        return SharedStateReference.untyped_reference(self._dictionary, full_key)

    @overload  # type: ignore
    def get(
        self, key: str, default: Optional[SharedStateReference[Serializable]] = None, /
    ) -> SharedStateReference[Serializable]:  # noqa: D102
        pass

    @overload
    def get(
        self,
        key: str,
        default: Optional[SharedStateReference[TSerializableObjectType]] = None,
        *,
        snapshot_type: Type[TSerializableObjectType],
        # noqa: D102
    ) -> SharedStateReference[TSerializableObjectType]:
        pass

    @override
    def get(  # type: ignore[override]
        self,
        key: str,
        default: Optional[SharedStateReference] = None,
        *,
        snapshot_type: Optional[Type[TSerializableObjectType]] = None,
    ) -> SharedStateReference:
        """Get a reference to a specific key."""
        if snapshot_type is None:
            return SharedStateReference.untyped_reference(self._dictionary, key)
        else:
            return SharedStateReference.typed_reference(
                self._dictionary, key, snapshot_type
            )
