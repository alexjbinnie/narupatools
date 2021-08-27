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

"""Utility methods for dealing with MDAnalysis."""

import contextlib

import numpy as np
from MDAnalysis import AtomGroup, NoDataError
from MDAnalysis.topology.guessers import guess_atom_element
from MDAnalysis.topology.tables import SYMB2Z, masses


def guess_elements(group: AtomGroup) -> np.ndarray:
    """
    Get the list of atomic symbols from a MDAnalysis object.

    If the group does not have element information, try converting atom names if
    present. If these are also not present, try using the atomic masses if available.

    :param group: AtomGroup to extract elements from.
    :raises NoDataError: AtomGroup has no elements, atom names or atom masses.
    """
    with contextlib.suppress(NoDataError):
        return group.elements
    with contextlib.suppress(NoDataError):
        return _atom_names_to_elements(group.names)
    with contextlib.suppress(NoDataError):
        return _atom_masses_to_elements(group.masses)
    raise NoDataError


def guess_atomic_number(group: AtomGroup) -> np.ndarray:
    """
    Get the list of atomic numbers from an MDAnalysis object.

    In order, this will try to convert elements, atom names and masses to a list of
    atomic numbers, skipping if that field is not available.

    :param group: AtomGroup to extract atomic numbers from.
    :raises NoDataError: AtomGroup has no elements, atom names or atom masses.
    """
    with contextlib.suppress(NoDataError):
        return _elements_to_atomic_number(group.elements)
    with contextlib.suppress(NoDataError):
        return _atom_names_to_atomic_number(group.names)
    with contextlib.suppress(NoDataError):
        return _atom_masses_to_atomic_number(group.masses)
    raise NoDataError


def _atom_mass_to_element(mass: float) -> str:
    for key, value in masses.items():
        if abs(value - mass) < 0.1:
            return key
    raise ValueError(f"Failed to guess Element from mass {mass}")


def _atom_names_to_atomic_number(names: np.ndarray) -> np.ndarray:
    return np.array(
        [SYMB2Z[guess_atom_element(name).capitalize()] for name in names], dtype=int
    )


def _elements_to_atomic_number(elements: np.ndarray) -> np.ndarray:
    return np.array([SYMB2Z[element] for element in elements], dtype=int)


def _atom_masses_to_atomic_number(elements: np.ndarray) -> np.ndarray:
    return np.array(
        [SYMB2Z[_atom_mass_to_element(element)] for element in elements], dtype=int
    )


def _atom_names_to_elements(names: np.ndarray) -> np.ndarray:
    return np.array(
        [guess_atom_element(name).capitalize() for name in names], dtype=object
    )


def _atom_masses_to_elements(elements: np.ndarray) -> np.ndarray:
    return np.array(
        [_atom_mass_to_element(element) for element in elements], dtype=object
    )
