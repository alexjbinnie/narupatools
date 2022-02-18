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

"""Dynamic NGLWidget that shows narupatools client."""

from typing import Any, Optional

from nglview import NGLWidget
from nglview.component import ComponentViewer

from narupatools.app import Client
from narupatools.frame.fields import ParticlePositions
from narupatools.nglview._structure import FrameDataStructure
from narupatools.util.timing import throttle


class _ClientWidget:
    def __init__(self, client: Client):
        self.widget = NGLWidget()
        self.client = client
        self.frame_component: Optional[ComponentViewer] = None
        if client.current_frame is not None:
            self._add_structure()
        client.on_frame_received.add_callback(self._on_frame)

    def _on_frame(self, **kwargs: Any) -> None:
        self.refresh()

    @throttle(0.05)
    def refresh(self) -> None:
        frame = self.client.current_frame
        if self.frame_component is None:
            self.frame_component = self._add_structure()
        positions = ParticlePositions.get(frame)
        self.frame_component.set_coordinates(10.0 * positions)

    def show(self) -> NGLWidget:
        return self.widget

    def _add_structure(self) -> ComponentViewer:
        frame = self.client.current_frame
        structure = FrameDataStructure(frame)
        return self.widget.add_structure(structure)


def show_client(client: Client) -> NGLWidget:
    """Show an NGLWidget that is tied to a Narupa client."""
    return _ClientWidget(client).show()
