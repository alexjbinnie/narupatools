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
from openmm import CustomIntegrator

class AMDIntegrator(CustomIntegrator):
    def __init__(self, dt: float, alpha, E) -> None: ...
    def getAlpha(self) -> None: ...
    def getE(self) -> None: ...
    def getEffectiveEnergy(self, energy) -> None: ...
    def setAlpha(self, alpha) -> None: ...
    def setE(self, E) -> None: ...

class AMDForceGroupIntegrator(CustomIntegrator):
    def __init__(self, dt: float, group, alphaGroup, EGroup) -> None: ...
    def getAlphaGroup(self) -> None: ...
    def getEGroup(self) -> None: ...
    def getEffectiveEnergy(self, groupEnergy: float) -> float: ...
    def setAlphaGroup(self, alpha) -> None: ...
    def setEGroup(self, E) -> None: ...

class DualAMDIntegrator(CustomIntegrator):
    def __init__(
        self, dt: float, group, alphaTotal, ETotal, alphaGroup, EGroup
    ) -> None: ...
    def getAlphaGroup(self) -> None: ...
    def getAlphaTotal(self) -> None: ...
    def getEGroup(self) -> None: ...
    def getETotal(self) -> None: ...
    def getEffectiveEnergy(self, totalEnergy: float, groupEnergy) -> float: ...
    def setAlphaGroup(self, alpha) -> None: ...
    def setAlphaTotal(self, alpha) -> None: ...
    def setEGroup(self, E) -> None: ...
    def setETotal(self, E) -> None: ...
