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
from ase import Atoms

from narupatools.ase.constraints import ASEConstraint


class NullConstraint(ASEConstraint):
    """ASE constraint that does nothing."""

    def adjust_positions(  # noqa: D102
        self, atoms: Atoms, positions: np.ndarray, /
    ) -> None:
        pass

    def adjust_forces(self, atoms: Atoms, forces: np.ndarray, /) -> None:  # noqa: D102
        pass
