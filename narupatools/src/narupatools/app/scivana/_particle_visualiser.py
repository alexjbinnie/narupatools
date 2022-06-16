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

from typing import Any, Union

import numpy as np
from narupa.trajectory import FrameData

import narupatools.util.properties as properties
from narupatools.state import SharedStateObject


class ParticleVisualisation(SharedStateObject):
    """Visualiation of a set of particles."""

    @properties.string
    def display_name(self) -> str:
        """User-facing display name."""

    @properties.auto
    def selection(self) -> Union[str, np.ndarray]:
        """
        Selection to be visualised, as either a reference to a predefined selection or a list of particle ids.

        If not set, all particles will be rendered.
        """

    @properties.auto
    def renderer(self) -> Any:
        """Renderer to be used to render the system."""

    @properties.boolean
    def hide(self) -> bool:
        """Should the visualiser be hidden."""

    @properties.integer
    def layer(self) -> int:
        """Layer on which the visualisation should be shown."""

    @properties.number
    def priority(self) -> float:
        """Priority of the visualisation in its layer."""

    @properties.auto
    def frame(self) -> FrameData:
        """
        Custom frame data to be used with the visualisaiton.

        If extend is set to True, then this frame data is overlaid over the existing frame data.

        If extend is set to False, then this frame data replaces the existing frame data. This is
        to allow ghost atoms to be shown.
        """

    @properties.boolean
    def extend(self) -> bool:
        """Should the additional frame data extend or replace the default?"""
