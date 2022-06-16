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
from typing import Any, Optional

import numpy as np
from ase.cell import Cell
from ase.geometry import find_mic
from narupa.trajectory import FrameData

from narupatools.frame import (
    BoxVectors,
    ParticleCount,
    ParticleElements,
    ParticlePositions,
    select,
)
from narupatools.physics.atomic import vdw_radius


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


def calculate_periodic_bonds(
    frame: FrameData,
    selection: str = None,
    method: Any = vdw_radius,
    cutoff_factor: float = 0.55,
) -> np.ndarray:
    """
    Calculate bonds in a periodic cell.

    :param frame: Molecular frame to use.
    :param selection: Selection string for which atoms to calculate bonds for.
    :param method: Method that returns a radius given an element.
    :param cutoff_factor: Cutoff factor. Bonds are valid if length < factor * (r(A) + r(B))
    :return: Bond pairs.
    """
    positions = frame[ParticlePositions]
    elements = frame[ParticleElements]
    cell = Cell(frame[BoxVectors])
    if selection is not None:
        selection_indices = select(frame, selection)
        positions = positions[selection_indices]
        elements = elements[selection_indices]
    indices = np.array(np.triu_indices(len(positions), k=1)).T
    offsets = positions[indices[:, 1]] - positions[indices[:, 0]]
    distances = find_mic(offsets, cell)[1]
    bond_cutoffs = cutoff_factor * method(elements[indices]).sum(axis=-1)
    if selection is not None:
        return selection_indices[indices[distances < bond_cutoffs]]  # type: ignore
    else:
        return indices[distances < bond_cutoffs]  # type: ignore


def supercell(
    frame: FrameData,
    n: Optional[int] = None,
    *,
    nx: Optional[int] = None,
    ny: Optional[int] = None,
    nz: Optional[int] = None
) -> FrameData:
    """Extend a frame to encompass a supercell."""
    if nx is None:
        nx = n
    if ny is None:
        ny = n
    if nz is None:
        nz = n

    if nx is None or ny is None or nz is None:
        raise ValueError("Either nx/ny/nz or n must be specified.")

    if nx < 1 or ny < 1 or nz < 1:
        raise ValueError("All dimensions must be >= 1.")

    total = nx * ny * nz - 1

    if ParticleCount in frame:
        frame[ParticleCount] *= total

    if BoxVectors not in frame:
        raise ValueError("Box vectors  missing.")

    vectors = frame[BoxVectors]
    x_axis = vectors[0]
    y_axis = vectors[1]
    z_axis = vectors[2]

    for key in frame:  # type: ignore
        if key.startswith("particle.") and key != "particle.count":
            if key == ParticlePositions:
                positions = frame[key]
                atom_count = len(positions)
                positions = np.concatenate([positions] * nx * ny * nz, axis=0)
                for i in range(nx):
                    for j in range(ny):
                        for k in range(nz):
                            print(
                                atom_count * (i + j * nx + k * (nx + ny)),
                                atom_count * (1 + i + j * nx + k * (nx + ny)),
                            )
                            positions[
                                atom_count
                                * (i + j * nx + k * (nx + ny)) : atom_count
                                * (1 + i + j * nx + k * (nx + ny)),
                                :,
                            ] += (
                                i * x_axis + j * y_axis + k * z_axis
                            )
                frame[key] = positions
            else:
                frame[key] = np.concatenate([frame[key]] * nx * ny * nz, axis=0)

    frame[BoxVectors] = np.array([x_axis * nx, y_axis * ny, z_axis * nz])

    return frame
