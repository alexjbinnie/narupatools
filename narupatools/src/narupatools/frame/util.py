"""Utility methods."""

import numpy as np


def calculate_residue_entities(
    *, residue_count: int, particle_residues: np.ndarray, bond_pairs: np.ndarray
) -> np.ndarray:
    """
    Calculate entities based on bonds between residues.

    An entity is a set of one or more residues which are connected by bonds.
    """
    entities = np.arange(residue_count)
    for pair in particle_residues[bond_pairs]:
        entities[pair.max()] = entities[pair.min()]

    return np.unique(entities, return_inverse=True)[1]  # type: ignore
