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

from os import PathLike
from typing import IO, Sequence, TextIO, Union, overload

import numpy as np
from simtk.openmm import Vec3
from simtk.openmm.app import Topology
from simtk.unit import Quantity
from typing_extensions import Literal

class PDBFile:
    topology: Topology
    positions: Sequence[Vec3]
    def __init__(self, *args: Union[bytes, str, PathLike, IO]): ...
    def getTopology(self) -> Topology: ...
    @overload
    def getPositions(self, asNumpy: Literal[True]) -> Quantity[np.ndarray]: ...
    @overload
    def getPositions(self, asNumpy: Literal[False]) -> Quantity[Sequence[Vec3]]: ...
    @overload
    def getPositions(
        self, asNumpy: bool
    ) -> Union[Quantity[np.ndarray], Quantity[Sequence[Vec3]]]: ...
    @staticmethod
    def writeFile(
        topology: Topology,
        positions: Union[Quantity[Sequence[Vec3]], Quantity[np.ndarray]],
        file: TextIO = ...,
    ) -> None: ...
