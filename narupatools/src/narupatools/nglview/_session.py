"""Dynamic NGLWidget that shows dynamics."""

from typing import Any

import numpy as np
from infinite_sets import everything, InfiniteSet
from nglview import NGLWidget

from narupatools.app import Session
from narupatools.app.scivana import CameraView
from narupatools.core.dynamics import SimulationDynamics
from narupatools.frame.fields import ParticlePositions, ParticleNames, BondPairs, ParticleResidues, ParticleElements, \
    ResidueNames, ResidueChains
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
    def __init__(self, session: Session):
        self.widget = NGLWidget()
        self.session = session
        frame = session.get_frame(fields=everything())
        structure = FrameDataStructure(frame)
        self.frame_component = self.widget.add_structure(structure)
        session.on_fields_changed.add_callback(self._on_fields_changed)
        self.widget.observe(self._on_camera_orientation_changed, names=["_camera_orientation"])
        self._camera_view = self.session.camera_views.add(CameraView(display_name="nglview"))
        self.update_camera(self.widget._camera_orientation)

    def _on_camera_orientation_changed(self, changes):
        self.update_camera(changes["new"])

    def update_camera(self, transformation):
        transformation = np.array(transformation)
        if len(transformation) != 16:
            return
        transformation = transformation[[0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14]]
        self._camera_view.update(transformation=transformation)

    def _on_fields_changed(self, *, fields: InfiniteSet[str] = everything(), **kwargs: Any) -> None:
        self.refresh(fields)

    @throttle(0.05)
    def refresh(self, fields: InfiniteSet[str]) -> None:
        changed_fields = _PDB_FIELDS & fields
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


def show_session(session: Session, /) -> NGLWidget:
    """
    Create an NGLWidget that dynamically shows a narupatools :class:`Session`.

    The widget automatically listens to changes in the underlying session and updates accordingly.
    """
    return _SessionWidget(session).show()
