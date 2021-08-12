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

"""Hard-coded globals defined by LAMMPS."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, TypeVar, Union, overload

from narupatools.lammps.constants import VariableType
from narupatools.lammps.wrapper import Extractable, LAMMPSWrapper

_TReturnType = TypeVar("_TReturnType")


@dataclass
class Global(Extractable[_TReturnType]):
    """A global value defined by LAMMPS that can be extracted using `extract_global`."""

    key: str

    def extract(self, lammps: LAMMPSWrapper) -> _TReturnType:
        """Extract the current value from a LAMMPS instance."""
        return lammps.extract_global(self.key)  # type: ignore[return-value]

    @overload
    @classmethod
    def define(cls, key: str, type: Literal[VariableType.INTEGER]) -> Global[int]:
        ...

    @overload
    @classmethod
    def define(cls, key: str, type: Literal[VariableType.DOUBLE]) -> Global[float]:
        ...

    @overload
    @classmethod
    def define(cls, key: str, type: Literal[VariableType.STRING]) -> Global[str]:
        ...

    @classmethod
    def define(
        cls,
        key: str,
        type: Literal[VariableType.DOUBLE, VariableType.INTEGER, VariableType.STRING],
    ) -> Union[Global[str], Global[int], Global[float]]:
        """Define a global that is defined by LAMMPS."""
        return Global(key)  # type: ignore[return-value]


TimestepLength = Global.define("dt", VariableType.DOUBLE)
""""Length of the timestep, in internal LAMMPS units."""

CurrentTimestep = Global.define("ntimestep", VariableType.INTEGER)
"""Current timestep as stored by LAMMPS."""

AccumulatedSimulationTime = Global.define("atime", VariableType.DOUBLE)
"""Accumulated simulation time in internal LAMMPS units."""

LastAccumulatedTimestep = Global.define("atimestep", VariableType.INTEGER)
"""The number of the timestep when the accumulated simulation time was last updated."""

UnitStyle = Global.define("units", VariableType.STRING)
"""The unit style currently in use."""
