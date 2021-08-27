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

"""Constants used by the LAMMPS API."""
from ctypes import c_double, c_int
from enum import IntEnum
from typing import Any

from lammps import LAMMPS_DOUBLE, LAMMPS_INT


class VariableStyle(IntEnum):
    """Variable styles that can be used for extract_compute."""

    GLOBAL = 0
    """Global variable."""
    ATOM = 1
    """Per atom variable."""
    LOCAL = 2
    """Local variable."""


class PropertyType(IntEnum):
    """Property styles that can be used for get_atom_property."""

    INT = LAMMPS_INT
    """Property consists of one or more integers."""
    DOUBLE = LAMMPS_DOUBLE
    """Property consists of one or more floating point numbers."""

    def get_ctype(self) -> Any:
        """Get the C type represented by this property type."""
        if self == PropertyType.INT:
            return c_int
        elif self == PropertyType.DOUBLE:
            return c_double


class VariableType(IntEnum):
    """Variable types that can be used for extract_compute."""

    SCALAR = 0
    """Variable is a single scalar value."""
    VECTOR = 1
    """Variable is a 1D array."""
    ARRAY = 2
    """Variable is a 2D array."""
