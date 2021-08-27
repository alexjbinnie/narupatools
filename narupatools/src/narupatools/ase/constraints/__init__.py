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
Code relating to ASE constraints, which can modify an ASE atoms object.

ASE constraints don't exist as a base class from which all other constraints subclass.
Instead, it is assumed that they have at a minimum an `adjust_positions` and
`adjust_forces` method, with optional other methods depending on what they alter.

This module includes protocols that provide these methods as abstract methods, ensuring
that if you include them in your class hierarchy then you have to implement them correctly.
"""

from ._constraint import (
    ASECellConstraint,
    ASEConstraint,
    ASEEnergyConstraint,
    ASEMomentaConstraint,
    ASEStressConstraint,
)
from ._imd_constraint import InteractionConstraint
from ._observer import ASEObserver

__all__ = [
    "InteractionConstraint",
    "ASEObserver",
    "ASEConstraint",
    "ASEMomentaConstraint",
    "ASEEnergyConstraint",
    "ASECellConstraint",
    "ASEStressConstraint",
]
