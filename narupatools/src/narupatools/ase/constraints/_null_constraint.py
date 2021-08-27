import numpy as np
from ase import Atoms

from narupatools.ase.constraints import ASEConstraint


class NullConstraint(ASEConstraint):
    def adjust_positions(self, atoms: Atoms, positions: np.ndarray, /) -> None:
        pass

    def adjust_forces(self, atoms: Atoms, forces: np.ndarray, /) -> None:
        pass