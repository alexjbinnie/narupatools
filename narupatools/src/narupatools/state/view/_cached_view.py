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

from typing import AbstractSet, Any, Dict, Generic, Optional, Set, Type, TypeVar, Union

from narupatools.override import override

from narupatools.state.typing import Serializable, SerializableObject
from ._collection_view import SharedStateCollectionView
from ._tracked_view import TrackedSharedStateView

TValue = TypeVar("TValue", bound=Union[SerializableObject, Serializable])


class CachedSharedStateCollectionView(
    SharedStateCollectionView[TValue], Generic[TValue]
):
    """View of a collection in the shared state, but caches the values extracted from it."""

    def __init__(
        self,
        session: TrackedSharedStateView,
        prefix: str,
        snapshot_type: Optional[Type[SerializableObject]] = None,
    ):
        self._values: Dict[str, TValue] = {}
        session.on_shared_state_added.add_callback(self._on_key_changed)
        session.on_shared_state_changed.add_callback(self._on_key_changed)
        session.on_shared_state_removed.add_callback(self._on_key_removed)
        self._changed_keys: Set[str] = set()
        self._removed_keys: Set[str] = set()
        super().__init__(session._dictionary, prefix, snapshot_type)

    def _on_key_changed(self, key: str, **kwargs: Any):
        if key.startswith(self.prefix):
            self._changed_keys.add(key)

    def _on_key_removed(self, key: str, **kwargs: Any):
        if key.startswith(self.prefix):
            self._removed_keys.add(key)

    def _ensure_updated(self):
        for removed in self._removed_keys.copy():
            del self._values[removed]
        self._removed_keys.clear()
        for changed in self._changed_keys.copy():
            self._values[changed] = self._snapshot_class.deserialize(
                self._dictionary[changed]
            )
        self._changed_keys.clear()

    @override(SharedStateCollectionView.snapshot)
    def snapshot(self):  # noqa: D102
        self._ensure_updated()
        return self._values

    @override(SharedStateCollectionView._keys)
    def _keys(self) -> AbstractSet[str]:
        self._ensure_updated()
        return self._values.keys()

    def __len__(self):
        self._ensure_updated()
        return len(self._values)
