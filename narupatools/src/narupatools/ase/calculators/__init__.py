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

"""
Various ASE calculators for debugging and simple systems.

ASE delegates the calculation of energy, forces, charges and few other properties to a
'calculator', which can be equated to a forcefield in other engines.
"""

from .constant_calculator import ConstantCalculator
from .null_calculator import NullCalculator
from .onebody_potential_calculator import OneBodyPotentialCalculator
from .protocols import CalculatorSetAtoms

__all__ = [
    "ConstantCalculator",
    "OneBodyPotentialCalculator",
    "CalculatorSetAtoms",
    "NullCalculator",
]
