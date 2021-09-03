"""Dynamic NGLWidget that shows narupatools client."""

from typing import Any

from nglview import NGLWidget

from narupatools.app import Client
from narupatools.core.timing import throttle
from narupatools.frame.fields import ParticlePositions
from narupatools.nglview._structure import FrameDataStructure


class _ClientWidget:
    def __init__(self, client: Client):
        self.widget = NGLWidget()
        self.client = client
        self.frame_component = None
        if client.current_frame is not None:
            self._add_structure()
        client.on_frame_received.add_callback(self._on_frame)

    def _on_frame(self, **kwargs: Any) -> None:
        self.refresh()

    @throttle(0.05)
    def refresh(self) -> None:
        frame = self.client.current_frame
        if self.frame_component is None:
            self._add_structure()
        positions = ParticlePositions.get(frame)
        self.frame_component.set_coordinates(10.0 * positions)

    def show(self) -> NGLWidget:
        return self.widget

    def _add_structure(self):
        frame = self.client.current_frame
        structure = FrameDataStructure(frame)
        self.frame_component = self.widget.add_structure(structure)


def show_client(client: Client) -> NGLWidget:
    return _ClientWidget(client).show()
