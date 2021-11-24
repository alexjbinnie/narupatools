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

from abc import ABCMeta, abstractmethod
from typing import Any

import numpy as np
from ase import Atom, Atoms

from narupatools.ase.calculators._calculator import Calculator
from narupatools.override import override
from narupatools.physics.typing import Vector3


class OneBodyPotentialCalculator(Calculator, metaclass=ABCMeta):
    """ASE calculator whose potential depends on each atom individually."""

    implemented_properties = ["forces", "energy"]

    @override(Calculator._calculate)
    def _calculate(self, atoms: Atoms, **kwargs: Any) -> None:  # noqa: D102
        forces = np.zeros((len(atoms), 3))
        energy = 0.0
        for i, atom in enumerate(atoms):
            energy += self.calculate_energy(atom)
            forces[i] += self.calculate_force(atom)
        self.results = {"forces": forces, "energy": energy}

    @abstractmethod
    def calculate_energy(self, atom: Atom) -> float:
        """
        Calculate the potential energy for a given atom, in electronvolts.

        :param atom: Atom to calculate the energy for.
        :return: Potential energy in electronvolts.
        """
        raise NotImplementedError

    @abstractmethod
    def calculate_force(self, atom: Atom) -> Vector3:
        """
        Calculate the force on a given atom, in electronvolts per angstrom.

        :param atom: Atom to calculate the force for.
        :return: Force in electronvolts per angstrom.
        """
        raise NotImplementedError
