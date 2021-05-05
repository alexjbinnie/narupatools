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

"""ASE calculator that is initialized with constant values."""

from typing import Any, Collection, List, Optional

import numpy as np
from ase.atoms import Atoms
from ase.calculators.calculator import Calculator, CalculatorSetupError, all_changes

from narupatools.physics.typing import ScalarArrayLike, Vector3ArrayLike
from .protocols import CalculatorSetAtoms


class ConstantCalculator(Calculator, CalculatorSetAtoms):
    """
    ASE calculator which is initialized with constant values.

    The calculator can be provided with a set of forces, set of charges or a potential
    energy. The calculator will always return this value, regardless of the system.
    However, an exception will be raised if the calculator is used with a system where
    the provided forces or charges have a different length to the size of the system.
    """

    def __init__(
        self,
        *,
        forces: Optional[Vector3ArrayLike] = None,
        energy: Optional[float] = None,
        charges: Optional[ScalarArrayLike] = None,
        **kwargs: Any,
    ):
        """
        Create a new calculator that uses the constant values provided.

        :param forces: Forces in electronvolts per angstroms.
        :param energy: Energy in electronvolts.
        :param charges: Charges in elementary charges.
        :param kwargs: Keyword arguments for the base
                       :class:`~ase.calculators.calculator.Calculator`.
        """
        super().__init__(**kwargs)
        self.implemented_properties = []
        if forces is not None:
            self.results["forces"] = np.asfarray(forces)
            self.implemented_properties.append("forces")
        if energy is not None:
            self.results["energy"] = energy
            self.implemented_properties.append("energy")
        if charges is not None:
            self.results["charges"] = np.asfarray(charges)
            self.implemented_properties.append("charges")

        if self.atoms is not None:
            self._check_atoms_is_valid(self.atoms)

    def reset(self) -> None:  # noqa: D102
        pass

    def calculate(  # noqa: D102
        self,
        atoms: Optional[Atoms] = None,
        properties: Collection[str] = ("forces", "energy"),
        system_changes: List[str] = all_changes,
    ) -> None:
        if atoms is not None:
            self._check_atoms_is_valid(atoms)
        super().calculate(atoms, properties, system_changes)

    def set_atoms(self, atoms: Atoms) -> None:  # noqa: D102
        self._check_atoms_is_valid(atoms)
        self.atoms = atoms

    def _check_atoms_is_valid(self, atoms: Atoms) -> None:
        if "forces" in self.results and len(self.results["forces"]) != len(atoms):
            raise CalculatorSetupError(
                "`forces` array is not the same size as the system."
            )
        if "charges" in self.results and len(self.results["charges"]) != len(atoms):
            raise CalculatorSetupError(
                "`charges` array is not the same size as the system."
            )
