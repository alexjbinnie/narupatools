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

from __future__ import annotations

from typing import TypeVar, overload

import numpy as np
from openmm.unit.quantity import Quantity

_T = TypeVar("_T")

class Unit:
    def __truediv__(self, other: Unit) -> Unit: ...
    def __rtruediv__(self, other: _T) -> Quantity[_T]: ...
    def __pow__(self, power: int) -> Unit: ...
    @overload
    def __mul__(self, other: np.ndarray) -> Quantity[np.ndarray]: ...
    @overload
    def __mul__(self, other: float) -> Quantity[float]: ...
    @overload
    def __mul__(self, other: Unit) -> Unit: ...
    @overload
    def __rmul__(self, other: float) -> Quantity[float]: ...
    @overload
    def __rmul__(self, other: np.ndarray) -> Quantity[np.ndarray]: ...
    @overload
    def __rmul__(self, other: Unit) -> Unit: ...
