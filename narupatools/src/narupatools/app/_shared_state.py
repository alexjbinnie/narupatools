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

"""Shared state view for clients and servers."""

import time
from typing import Optional, Protocol, Set

from narupa.utilities.change_buffers import DictionaryChange

from narupatools.core.event import Event, EventListener
from narupatools.imd.interactions._interactiondata import InteractionData
from narupatools.state import SharedStateCollectionView, SharedStateDictionaryView
from narupatools.state.typing import SerializableDictionary


class OnSharedStateAddedCallback(Protocol):
    """Callback for when a key is first added to the shared state."""

    def __call__(self, *, key: str, token: Optional[str], time_added: float) -> None:
        """
        Called when a key is first added to the shared state.

        :param key: Full key of the item.
        :param token: Access token used to add this item.
        :param time_added: Real time that this key was added.
        """


class OnSharedStateChangedCallback(Protocol):
    """Callback for when a key is changed in the shared state."""

    def __call__(self, *, key: str, token: Optional[str], time_updated: float) -> None:
        """
        Called when a key is changed in the shared state.

        :param key: Full key of the item.
        :param token: Access token that is modifying this item.
        :param time_updated: Real time that this key was updated.
        """


class OnSharedStateRemovedCallback(Protocol):
    """Callback for when a key is removed from the shared state."""

    def __call__(self, *, key: str, token: Optional[str], time_deleted: float) -> None:
        """
        Called when a key is removed from the shared state.

        :param key: Full key of the item.
        :param token: Access token that is deleting this item.
        :param time_deleted: Real time that this key was deleted.
        """


class SessionSharedState(SharedStateDictionaryView):
    """View of a shared state with callbacks for when keys are changed."""

    def __init__(self, dictionary: SerializableDictionary):
        """
        Create a new view.

        :param dictionary: Shared state dictionary this is a view of.
        """
        super().__init__(dictionary)
        self._on_added = Event(OnSharedStateAddedCallback)
        self._on_changed = Event(OnSharedStateChangedCallback)
        self._on_removed = Event(OnSharedStateRemovedCallback)
        self._key_history: Set[str] = set()

    def on_shared_state_added(self) -> EventListener[OnSharedStateAddedCallback]:
        """Callback when a new item is added to the shared state."""
        return self._on_added

    def on_shared_state_changed(self) -> EventListener[OnSharedStateChangedCallback]:
        """Callback when an item is modified in the shared state."""
        return self._on_changed

    def on_shared_state_removed(self) -> EventListener[OnSharedStateRemovedCallback]:
        """Callback when an item is removed from the shared state."""
        return self._on_removed

    def _on_dictionary_update(
        self, *, access_token: Optional[str], change: DictionaryChange
    ) -> None:
        """
        Trigger callbacks based on changes to a shared state dictionary.

        :param access_token: Token used to modify the state.
        :param change: Changes that have occured.
        """
        current_time = time.monotonic()
        for key in change.updates:
            if key not in self._key_history:
                self._key_history.add(key)
                self._on_added.invoke(
                    key=key, token=access_token, time_added=current_time
                )
            else:
                self._on_changed.invoke(
                    key=key, token=access_token, time_updated=current_time
                )
        for key in change.removals:
            if key in self._key_history:
                self._key_history.remove(key)
                self._on_removed.invoke(
                    key=key, token=access_token, time_deleted=current_time
                )

    @property
    def interactions(self) -> SharedStateCollectionView:
        """View of current interactions affecting the system."""
        return self.collection("interaction", InteractionData)
