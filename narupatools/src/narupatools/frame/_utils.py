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
from typing import List, Optional

import numpy as np
from MDAnalysis.topology.tables import SYMB2Z, Z2SYMB, masses


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


__mass_cutoff_table: Optional[List[List[int]]] = None


def mass_to_element(mass: float, tolerance: float = 0.1) -> Optional[int]:
    """Get the atomic number with the closest mass to a given atomic mass."""
    global __mass_cutoff_table

    if __mass_cutoff_table is None:
        mass_dict = masses.copy()
        del [mass_dict["DUMMY"]]

        mass_list = np.array(list(mass_dict.values()))
        symbol_list = np.array(list(mass_dict.keys()))
        masses_sorted = np.argsort(mass_list)

        __mass_cutoff_table = []

        for i in range(len(masses_sorted) - 1):
            cutoff = 0.5 * (
                mass_list[masses_sorted[i]] + mass_list[masses_sorted[i + 1]]
            )
            __mass_cutoff_table.append(
                [
                    cutoff,
                    mass_list[masses_sorted[i]],
                    SYMB2Z[symbol_list[masses_sorted[i]].capitalize()],
                ]
            )

    for cutoff, actual, key in __mass_cutoff_table:
        if mass < cutoff:
            if abs(mass - actual) < tolerance:
                return key
            return None
    return None
