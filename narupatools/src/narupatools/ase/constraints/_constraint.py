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

from abc import abstractmethod
from typing import Protocol

import numpy as np
from ase import Atoms


class ASEConstraint(Protocol):
    """Protocol describing an ASE Constraint, as ASE uses duck typing."""

    @abstractmethod
    def adjust_positions(self, atoms: Atoms, positions: np.ndarray, /) -> None:
        """
        Adjust the positions in-place for an ASE `Atoms` object.

        :param atoms: The ASE `Atoms` object this constraint applies to.
        :param positions: The positions to be modified by this constraint, in angstrom.
        """
        raise NotImplementedError

    @abstractmethod
    def adjust_forces(self, atoms: Atoms, forces: np.ndarray, /) -> None:
        """
        Adjust the forces in-place for an ASE `Atoms` object.

        :param atoms: The ASE `Atoms` object this constraint applies to.
        :param forces: The forces to be modified by this constraint, in eV per nm.
        """
        raise NotImplementedError


class ASEEnergyConstraint(ASEConstraint, Protocol):
    """Protocol for a constraint that modifies potential energy."""

    @abstractmethod
    def adjust_potential_energy(self, atoms: Atoms, /) -> float:
        """
        Get the difference in potential energy due to this constraint.

        :param atoms: The ASE `Atoms` object this constraint applies to.
        :return: The difference in potential energy, in eV.
        """
        raise NotImplementedError


class ASEMomentaConstraint(ASEConstraint, Protocol):
    """Protocol for a Constraint that modifies the momenta."""

    @abstractmethod
    def adjust_momenta(self, atoms: Atoms, momenta: np.ndarray, /) -> None:
        """
        Adjust the momenta in-place for an ASE `Atoms` object.

        :param atoms: The ASE `Atoms` object this constraint applies to.
        :param momenta: The momenta to be modified by this constraint, in a.m.u angstrom
                        per ASE time unit.
        """
        raise NotImplementedError


class ASETorqueConstraint(ASEConstraint, Protocol):
    """Protocol for a Constraint that modifies the torques."""

    @abstractmethod
    def adjust_torques(self, atoms: Atoms, torques: np.ndarray, /) -> None:
        """
        Adjust the torques in-place for an ASE `Atoms` object.

        :param atoms: The ASE `Atoms` object this constraint applies to.
        :param torques: The torques to be modified by this constraint, in eV.
        """
        raise NotImplementedError


class ASECellConstraint(ASEConstraint, Protocol):
    """Protocol for a constraint that modifies the unit cell."""

    @abstractmethod
    def adjust_cell(self, atoms: Atoms, cell: np.ndarray, /) -> None:
        """
        Adjust the unit cell in-place for an ASE `Atoms` object.

        :param atoms:  The ASE `Atoms` object this constraint applies to.
        :param cell: 3x3 NumPy array of cell vectors, in angstrom.
        """
        raise NotImplementedError


class ASEStressConstraint(ASEConstraint, Protocol):
    """Protocol for a constraint that modifies stresses."""

    @abstractmethod
    def adjust_stress(self, atoms: Atoms, stress: np.ndarray, /) -> None:
        """
        Adjust the stress in-place for an ASE `Atoms` object.

        :param atoms:  The ASE `Atoms` object this constraint applies to.
        :param stress: Either a 6x1 (in Voigt order) or 3x3 NumPy array representing the
                       symmetric stress tensor, in electronvolts per nanometers cubed.
        """
        raise NotImplementedError
