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

from typing import Any, IO, Mapping, Optional, Union

from simtk.openmm import Context, Integrator, System
from simtk.openmm.app import Topology
from simtk.unit import Quantity

class Simulation:
    topology: Topology
    system: System
    integrator: Integrator
    currentStep: int
    context: Context
    def __init__(
        self,
        topology: Topology,
        system: Union[System, str],
        integrator: Union[Integrator, str],
        platform: Any = ...,
        platformProperties: Optional[Mapping[str, Any]] = ...,
        state: Optional[str] = ...,
    ): ...
    def minimizeEnergy(
        self, tolerance: Union[float, Quantity[float]] = ..., maxIterations: int = ...
    ) -> None: ...
    def step(self, steps: int) -> None: ...
    def runForClockTime(
        self,
        time: Union[float, Quantity[float]],
        checkpointFile: Optional[Union[str, IO]] = ...,
        stateFile: Optional[Union[str, IO]] = ...,
        checkpointInterval: Optional[Union[float, Quantity[float]]] = ...,
    ) -> None: ...
    def saveCheckpoint(self, file: Union[str, IO]) -> None: ...
    def loadCheckpoint(self, file: Union[str, IO]) -> None: ...
    def saveState(self, file: Union[str, IO]) -> None: ...
    def loadState(self, file: Union[str, IO]) -> None: ...
