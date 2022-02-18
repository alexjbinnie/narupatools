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
    import lammps

    from ._converter import LAMMPSConverter  # noqa: F401
    from ._dynamics import LAMMPSDynamics
    from ._simulation import LAMMPSSimulation

    lmp = lammps.lammps()
    INSTALLED_PACKAGES = set(lmp.installed_packages)
    """Set of packages currently installed with LAMMPS."""
    INSTALLED_ATOM_STYLES = set(lmp.available_styles("atom"))
    INSTALLED_INTEGRATE_STYLES = set(lmp.available_styles("integrate"))
    INSTALLED_MINIMIZE_STYLES = set(lmp.available_styles("minimize"))
    INSTALLED_PAIR_STYLES = set(lmp.available_styles("pair"))
    INSTALLED_BOND_STYLES = set(lmp.available_styles("bond"))
    INSTALLED_ANGLE_STYLES = set(lmp.available_styles("angle"))
    INSTALLED_DIHEDRAL_STYLES = set(lmp.available_styles("dihedral"))
    INSTALLED_IMPROPER_STYLES = set(lmp.available_styles("improper"))
    INSTALLED_KSPACE_STYLES = set(lmp.available_styles("kspace"))
    INSTALLED_COMPUTE_STYLES = set(lmp.available_styles("compute"))
    INSTALLED_FIX_STYLES = set(lmp.available_styles("fix"))
    INSTALLED_REGION_STYLES = set(lmp.available_styles("region"))
    INSTALLED_DUMP_STYLES = set(lmp.available_styles("dump"))
    INSTALLED_COMMAND_STYLES = set(lmp.available_styles("command"))

    lmp.close()


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
    "INSTALLED_PACKAGES",
    "INSTALLED_ATOM_STYLES",
    "INSTALLED_INTEGRATE_STYLES",
    "INSTALLED_MINIMIZE_STYLES",
    "INSTALLED_PAIR_STYLES",
    "INSTALLED_BOND_STYLES",
    "INSTALLED_ANGLE_STYLES",
    "INSTALLED_DIHEDRAL_STYLES",
    "INSTALLED_IMPROPER_STYLES",
    "INSTALLED_KSPACE_STYLES",
    "INSTALLED_COMPUTE_STYLES",
    "INSTALLED_FIX_STYLES",
    "INSTALLED_REGION_STYLES",
    "INSTALLED_DUMP_STYLES",
    "INSTALLED_COMMAND_STYLES",
]
