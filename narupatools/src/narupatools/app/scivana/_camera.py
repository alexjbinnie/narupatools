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

import numpy as np

from narupatools.state import SharedStateObject
from narupatools.util import properties


class CameraView(SharedStateObject):
    """Fixed camera view from a 2D application."""

    @properties.auto
    def transformation(self) -> np.ndarray:
        """Transformation of the camera as an affine transformation."""

    @properties.auto
    def display_name(self) -> str:
        """Display name of the camera."""
