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
from abc import abstractmethod
from typing import Optional, Union

import numpy as np
import numpy.typing as npt
from infinite_sets import everything
from MDAnalysis import Universe

from narupatools.app.scivana import (
    CameraView,
    ParticleSelection,
    ParticleVisualisation,
    Renderer,
)
from narupatools.frame import FrameSource, convert
from narupatools.imd.interactions import InteractionParameters
from narupatools.state import SharedStateReference
from narupatools.state.view import (
    CachedSharedStateCollectionView,
    TrackedSharedStateView,
)


class SharedStateMixin(FrameSource):
    """Mixin that adds common shared state operations to clients and sessions."""

    def __init__(self):
        self.interactions = CachedSharedStateCollectionView(
            self.shared_state, "interaction.", InteractionParameters
        )
        self.selections = CachedSharedStateCollectionView(
            self.shared_state, "selection.", ParticleSelection
        )
        self.visualisations = CachedSharedStateCollectionView(
            self.shared_state, "visualisation.", ParticleVisualisation
        )
        self.camera_views = CachedSharedStateCollectionView(
            self.shared_state, "camera.", CameraView
        )

    @property
    @abstractmethod
    def shared_state(self) -> TrackedSharedStateView:
        """View of the shared state."""

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
    ) -> SharedStateReference[ParticleVisualisation]:
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
        vis = ParticleVisualisation(
            selection=selection, renderer=renderer, layer=layer, priority=priority
        )
        return self.visualisations.add(vis)
