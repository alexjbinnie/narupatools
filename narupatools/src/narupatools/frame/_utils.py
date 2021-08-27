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

"""Utilities for converting different arrays."""

import numpy as np
from MDAnalysis.topology.tables import Z2SYMB, masses


def atomic_numbers_to_symbols(numbers: np.ndarray, /) -> np.ndarray:
    """
    Convert an array of atomic numbers to an array of element symbols.

    :param numbers: Array of atomic numbers to convert.
    """
    return np.array([Z2SYMB[number] for number in numbers])


def atomic_numbers_to_masses(numbers: np.ndarray, /) -> np.ndarray:
    """
    Convert an array of atomic numbers to an array of atomic masses.

    :param numbers: Array of atomic numbers to convert.
    """
    return np.array([masses[Z2SYMB[number]] for number in numbers])
