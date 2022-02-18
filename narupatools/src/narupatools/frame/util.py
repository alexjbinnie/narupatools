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

"""Utility methods."""

import numpy as np


def calculate_residue_entities(
    *, residue_count: int, particle_residues: np.ndarray, bond_pairs: np.ndarray
) -> np.ndarray:
    """
    Calculate entities based on bonds between residues.

    An entity is a set of one or more residues which are connected by bonds.
    """
    assert len(bond_pairs) > 0
    entities = np.arange(residue_count)
    for pair in particle_residues[bond_pairs]:
        entities[pair.max()] = entities[pair.min()]

    return np.unique(entities, return_inverse=True)[1]  # type: ignore
