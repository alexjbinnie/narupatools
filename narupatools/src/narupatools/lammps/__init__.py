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

"""Interpolability with the LAMMPS package."""

import importlib

from .exceptions import LAMMPSError, LAMMPSWarning

has_lammps = importlib.util.find_spec("lammps") is not None

if has_lammps:
    from ._converter import LAMMPSConverter  # noqa: F401
    from ._dynamics import LAMMPSDynamics
    from ._simulation import LAMMPSSimulation

from ._datafile import LAMMPSDataFile
from ._units import (
    UnitsLAMMPSCGS,
    UnitsLAMMPSElectron,
    UnitsLAMMPSMetal,
    UnitsLAMMPSMicro,
    UnitsLAMMPSNano,
    UnitsLAMMPSReal,
    UnitsLAMMPSSI,
    get_unit_system,
)

__all__ = [
    "UnitsLAMMPSCGS",
    "UnitsLAMMPSSI",
    "UnitsLAMMPSMicro",
    "UnitsLAMMPSMetal",
    "UnitsLAMMPSElectron",
    "UnitsLAMMPSNano",
    "UnitsLAMMPSReal",
    "get_unit_system",
    "LAMMPSSimulation",
    "LAMMPSError",
    "LAMMPSWarning",
    "LAMMPSDynamics",
    "LAMMPSDataFile",
]
