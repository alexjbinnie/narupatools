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

if not has_lammps:
    raise ImportError("narupatools.lammps requires lammps to be installed.")

from .calculator import LAMMPSCalculator
from .converter import atoms_from_lammps_simulation
from .dynamics import LAMMPSDynamics
from .simulation import LAMMPSSimulation
from .units import (
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
    "LAMMPSCalculator",
    "LAMMPSSimulation",
    "LAMMPSError",
    "LAMMPSWarning",
    "LAMMPSDynamics",
    "atoms_from_lammps_simulation",
]
