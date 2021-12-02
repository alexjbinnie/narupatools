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

from infinite_sets import everything
from nglview import NGLWidget

from narupatools.core.dynamics import SimulationDynamics
from narupatools.frame.fields import ParticlePositions
from narupatools.nglview._structure import FrameDataStructure
from narupatools.util.timing import throttle


class _DynamicsWidget:
    def __init__(self, dynamics: SimulationDynamics):
        self.widget = NGLWidget()
        self.dynamics = dynamics
        frame = dynamics.get_frame(fields=everything())
        structure = FrameDataStructure(frame)
        self.frame_component = self.widget.add_structure(structure)
        dynamics.on_fields_changed.add_callback(self._on_fields_changed)

    def _on_fields_changed(self, **kwargs: Any) -> None:
        self.refresh()

    @throttle(0.05)
    def refresh(self) -> None:
        frame = self.dynamics.get_frame(fields={ParticlePositions.key})
        positions = ParticlePositions.get(frame)
        self.frame_component.set_coordinates(10.0 * positions)

    def show(self) -> NGLWidget:
        return self.widget


def show_dynamics(dynamics: SimulationDynamics) -> NGLWidget:
    """
    Create an NGLWidget that dynamically shows a narupatools :class:`SimulationDynamics`.

    The widget automatically listens to changes in the underlying dynamics and updates accordingly.

    :param dynamics: Dynamics to listen to.
    :return: NGLWidget dynamically showing dynamics.
    """
    return _DynamicsWidget(dynamics).show()
