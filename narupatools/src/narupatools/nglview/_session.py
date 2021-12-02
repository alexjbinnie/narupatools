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

"""Dynamic NGLWidget that shows dynamics."""

from typing import Any

import numpy as np
from infinite_sets import InfiniteSet, everything
from nglview import NGLWidget

from narupatools.app import Session
from narupatools.app.scivana import CameraView
from narupatools.frame.fields import (
    BondPairs,
    ParticleElements,
    ParticleNames,
    ParticlePositions,
    ParticleResidues,
    ResidueChains,
    ResidueNames,
)
from narupatools.nglview._structure import FrameDataStructure
from narupatools.util.timing import throttle

# Set of fields which are required to define a structure
_PDB_FIELDS = {
    ParticlePositions.key,
    ParticleNames.key,
    BondPairs.key,
    ParticleResidues.key,
    ParticleElements.key,
    ResidueNames.key,
    ResidueChains.key,
}


class _SessionWidget:
    def __init__(self, session: Session, sync_camera: bool = False):
        self.widget = NGLWidget()
        self.session = session
        frame = session.get_frame(fields=everything())
        structure = FrameDataStructure(frame)
        self.frame_component = self.widget.add_structure(structure)
        session.on_fields_changed.add_callback(self._on_fields_changed)
        if sync_camera:
            self.widget.observe(
                self._on_camera_orientation_changed, names=["_camera_orientation"]
            )
            self._camera_view = self.session.camera_views.add(
                CameraView(display_name="nglview")
            )
            self.update_camera(self.widget._camera_orientation)

    def _on_camera_orientation_changed(self, changes: Any) -> None:
        self.update_camera(changes["new"])

    def update_camera(self, transformation: Any) -> None:
        transformation = np.array(transformation)
        if len(transformation) != 16:
            return
        transformation = transformation[[0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14]]
        self._camera_view.update(transformation=transformation)

    def _on_fields_changed(
        self, *, fields: InfiniteSet[str] = everything(), **kwargs: Any
    ) -> None:
        self.refresh(fields)

    @throttle(0.05)
    def refresh(self, fields: InfiniteSet[str]) -> None:
        changed_fields = fields & _PDB_FIELDS
        if changed_fields == {ParticlePositions.key}:
            frame = self.session.get_frame(fields={ParticlePositions.key})
            positions = ParticlePositions.get(frame)
            self.frame_component.set_coordinates(10.0 * positions)
        else:
            frame = self.session.get_frame(fields=everything())
            structure = FrameDataStructure(frame)
            self.widget.remove_component(self.frame_component)
            self.frame_component = self.widget.add_structure(structure)

    def show(self) -> NGLWidget:
        return self.widget


def show_session(session: Session, /, *, sync_camera: bool = False) -> NGLWidget:
    """
    Create an NGLWidget that dynamically shows a narupatools :class:`Session`.

    The widget automatically listens to changes in the underlying session and updates accordingly.

    :param sync_camera: Should the NGL camera view be synced through the session. This can be used to use the NGL
                        widget to align screenshots to be taken from the VR client.
    """
    return _SessionWidget(session, sync_camera=sync_camera).show()
