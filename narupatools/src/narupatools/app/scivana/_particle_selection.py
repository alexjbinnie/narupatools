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

import narupatools.util.properties as properties
from narupatools.state import SharedStateObject


class ParticleSelection(SharedStateObject):
    """Representation of a named predefined selection."""

    @properties.numpy(dtype=int, shape=(None,))
    def particle_ids(self) -> np.ndarray:
        """List of particle ids included in the selection."""
        ...

    @properties.string
    def display_name(self) -> str:
        """User-facing display name of the selection."""
        ...
