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
from abc import abstractmethod
from typing import Optional, Protocol, Set, Union

import numpy as np
import numpy.typing as npt
from infinite_sets import everything
from MDAnalysis import Universe
from narupa.utilities.change_buffers import DictionaryChange

from narupatools.app.scivana import ParticleSelection, ParticleVisualiser, Renderer
from narupatools.core.event import Event, EventListener
from narupatools.frame import FrameSource, convert
from narupatools.imd.interactions import InteractionParameters
from narupatools.state import (
    SharedStateCollectionView,
    SharedStateDictionaryView,
    SharedStateReference,
)
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


class SharedStateMixin(FrameSource, Protocol):
    """Mixin that adds common shared state operations to clients and sessions."""

    @property
    @abstractmethod
    def shared_state(self) -> SessionSharedState:
        """View of the shared state."""
        pass

    def create_selection(
        self,
        particles: Union[str, npt.ArrayLike],
        /,
        *,
        display_name: str = "selection",
        key: Optional[str] = None,
    ) -> SharedStateReference[ParticleSelection]:
        """
        Create a new selection.

        :param particles: MDAnalysis selection string or array of particle indices.
        :param display_name: Display name to be shown in UI.
        :param key: Key to store the selection. If not provided, one will be generated.
        :returns: Reference to the selection just defined.
        """
        if isinstance(particles, str):
            frame = self.get_frame(fields=everything())
            universe = convert(frame, Universe)
            atom_group = universe.select_atoms(particles)
            particles = atom_group.indices
        particles = np.asarray(particles, dtype=int)
        selection = ParticleSelection(particles=particles, display_name=display_name)
        if key is None:
            return self.selections.add(selection)
        else:
            return self.selections.set(key, selection)

    def create_visualisation(
        self,
        *,
        selection: Optional[
            Union[str, npt.ArrayLike, SharedStateReference[ParticleSelection]]
        ] = None,
        renderer: Renderer,
        layer: Optional[int] = None,
        priority: Optional[int] = None,
    ) -> SharedStateReference[ParticleVisualiser]:
        """
        Create a new visualisation to draw atoms.

        :param selection: MDAnalysis selection string, array of particle indices or a reference to a
                          previously defined selection.
        :param renderer: Renderer to use.
        :param layer: Layer to add the renderer to. Renderers on the same layer occlude each other.
        :param priority: Priority of the renderer in the layer. Renderers with higher priority
        """
        if isinstance(selection, str):
            frame = self.get_frame(fields=everything())
            universe = convert(frame, Universe)
            atom_group = universe.select_atoms(selection)
            selection = np.asarray(atom_group.indices, dtype=int)
        elif isinstance(selection, SharedStateReference):
            selection = selection.key
        elif selection is None:
            selection = None
        else:
            selection = np.asarray(selection, dtype=int)
        vis = ParticleVisualiser(
            selection=selection, renderer=renderer, layer=layer, priority=priority
        )
        return self.visualisers.add(vis)

    @property
    def interactions(self) -> SharedStateCollectionView[InteractionParameters]:
        """View of current interactions affecting the system."""
        return self.shared_state.collection("interaction.", InteractionParameters)

    @property
    def selections(self) -> SharedStateCollectionView[ParticleSelection]:
        """View of current selections defined for the system."""
        return self.shared_state.collection("selection.", ParticleSelection)

    @property
    def visualisers(self) -> SharedStateCollectionView[ParticleVisualiser]:
        return self.shared_state.collection("visualisation.", ParticleVisualiser)
