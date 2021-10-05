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

from typing import Union

import numpy as np
from narupa.trajectory import FrameData

import narupatools.util.properties as properties
from narupatools.app.scivana import Renderer
from narupatools.state import SharedStateObject


class ParticleVisualiser(SharedStateObject):
    @properties.string
    def display_name(self) -> str:
        pass

    @properties.string
    def selection(self) -> Union[str, np.ndarray]:
        pass

    @properties.auto
    def renderer(self) -> Renderer:
        pass

    @properties.boolean
    def hide(self) -> bool:
        pass

    @properties.integer
    def layer(self) -> int:
        pass

    @properties.number
    def priority(self) -> float:
        pass

    @properties.auto
    def frame(self) -> FrameData:
        pass

    @properties.boolean
    def extend(self) -> bool:
        pass
