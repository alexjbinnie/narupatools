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

from enum import IntEnum

from lammps import (
    LAMMPS_DOUBLE,
    LAMMPS_DOUBLE_2D,
    LAMMPS_INT,
    LAMMPS_INT_2D,
    LAMMPS_STRING,
)


class VariableStyle(IntEnum):
    """Variable styles that can be used for extract_compute."""

    GLOBAL = 0
    """Global variable."""
    ATOM = 1
    """Per atom variable."""
    LOCAL = 2
    """Local variable."""


class VariableType(IntEnum):
    """Variable types that can be used for extracting data from LAMMPS."""

    INTEGER = LAMMPS_INT
    """LAMMPS 32-bit integer."""

    INTEGER_ARRAY = LAMMPS_INT_2D
    """LAMMPS 32-bit integer."""

    DOUBLE = LAMMPS_DOUBLE
    """LAMMPS 64-bit double."""

    DOUBLE_ARRAY = LAMMPS_DOUBLE_2D
    """LAMMPS 64-bit double."""

    STRING = LAMMPS_STRING
    """LAMMPS String."""


class VariableDimension(IntEnum):
    """Variable types that can be used for extract_compute."""

    SCALAR = 0
    """Variable is a single scalar value."""
    VECTOR1D = 1
    """Variable is a 1D array."""
    ARRAY2D = 2
    """Variable is a 2D array."""
