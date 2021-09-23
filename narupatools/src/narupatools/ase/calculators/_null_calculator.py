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

from typing import Any

import numpy as np
from ase import Atoms
from ase.calculators.calculator import all_changes

from ...override import override
from ._calculator import Calculator


class NullCalculator(Calculator):
    """Empty ASE calculator which generates zero values for forces and energies."""

    ignored_changes = all_changes

    implemented_properties = ["forces", "energy"]

    @override
    def _calculate(self, atoms: Atoms, **kwargs: Any) -> None:  # noqa: D102
        self.results = {"forces": np.zeros((len(atoms), 3)), "energy": 0.0}
