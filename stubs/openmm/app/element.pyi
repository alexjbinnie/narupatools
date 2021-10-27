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

from simtk.unit import Quantity

class Element:
    mass: Quantity[float]
    symbol: str
    atomic_number: int
    name: str
    @classmethod
    def getBySymbol(cls, name: str) -> Element: ...
    @classmethod
    def getByAtomicNumber(cls, atomic_number: int) -> Element: ...
    @classmethod
    def getByMass(cls, mass: Union[float, Quantity[float]]) -> Element: ...
