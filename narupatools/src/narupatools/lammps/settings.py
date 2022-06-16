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

"""Hard-coded global settings defined by LAMMPS."""


from dataclasses import dataclass

from narupatools.lammps._wrapper import Extractable, LAMMPSWrapper


@dataclass
class Setting(Extractable[int]):
    """Setting that can be obtained from LAMMPS using `extract_setting`."""

    key: str

    def extract(self, lammps: LAMMPSWrapper) -> int:
        """Extract the current value of this setting from a LAMMPS instance."""
        return lammps.extract_setting(self.key)


@dataclass
class BooleanSetting(Extractable[bool]):
    """Boolean setting that can be obtained from LAMMPS using `extract_setting`."""

    key: str

    def extract(self, lammps: LAMMPSWrapper) -> bool:
        """Extract the current value of this setting from a LAMMPS instance."""
        return lammps.extract_setting(self.key) > 0


SizeBigInt = Setting("bigint")
"""Size of the bigint integer type in bytes."""

SizeTagInt = Setting("tagint")
"""Size of the tagint integer type in bytes"""

SizeImageInt = Setting("imageint")
"""Size of the imageint integer type in bytes."""

Dimension = Setting("dimension")
"""Dimension of the system."""

IsBoxDefined = BooleanSetting("box_exist")
"""Is the simulation box defined?"""

IsBoxTriclinic = BooleanSetting("triclinic")
"""Is the simulation box triclinic?"""

NumberLocalAtoms = Setting("nlocal")
"""The number of 'owned' atoms of the current MPI rank."""

NumberGhostAtoms = Setting("nghost")
"""The number of 'ghost' atoms of the current MPI rank."""

NumberAllAtoms = Setting("nall")
"""The number of 'ghost' and 'owned' atoms of the current MPI rank."""

NumberAtomTypes = Setting("ntypes")
"""The number of atom types."""

NumberBondTypes = Setting("nbondtypes")
"""The number of bond types."""

NumberAngleTypes = Setting("nangletypes")
"""The number of angle types."""

NumberDihedralTypes = Setting("ndihedraltypes")
"""The number of dihedral types."""

NumberImproperTypes = Setting("nimpropertypes")
"""The number of improper types."""

AtomStylesIncludesMolecularTopology = BooleanSetting("molecule_flag")
"""Does the atom style includes molecular topology data?"""

AtomStylesIncludesCharges = BooleanSetting("q_flag")
"""Does the atom style includes point charges?"""

AtomStylesIncludesDipoles = BooleanSetting("mu_flag")
"""Does the atom style includes point dipoles?"""

AtomStylesIncludesPerAtomMasses = BooleanSetting("rmass_flag")
"""Does the atom style includes per-atom masses?"""

AtomStylesIncludesPerAtomRadii = BooleanSetting("radius_flag")
"""Does the atom style includes a per-atom radius?"""

AtomStylesCanRotate = BooleanSetting("sphere_flag")
"""Does the atom style describes extended particles that can rotate?"""

AtomStylesCanBeEllipsoid = BooleanSetting("ellipsoid_flag")
"""Does the atom style describes extended particles that may be ellipsoidal?"""

AtomStylesIncludesRotationalVelocities = BooleanSetting("omega_flag")
"""Does the atom style store per-atom rotational velocities?"""

AtomStylesIncludesTorques = BooleanSetting("torque_flag")
"""Does the atom style include per-atom torques?"""

AtomStylesIncludesAngularMomenta = BooleanSetting("angmom_flag")
"""Does the atom style can store per-atom angular momentum?"""
