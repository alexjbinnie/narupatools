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

from typing import Any, Optional

import numpy as np
from ase.atoms import Atoms
from ase.calculators.calculator import CalculatorSetupError, all_changes

from narupatools.ase.calculators._calculator import Calculator
from narupatools.physics.typing import ScalarArrayLike, Vector3ArrayLike


class ConstantCalculator(Calculator):
    """
    ASE calculator which is initialized with constant values.

    The calculator can be provided with a set of forces, set of charges or a potential
    energy. The calculator will always return this value, regardless of the system.
    However, an exception will be raised if the calculator is used with a system where
    the provided forces or charges have a different length to the size of the system.
    """

    ignored_changes = all_changes

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
            self._forces = np.asfarray(forces)
            self.implemented_properties.append("forces")
        if energy is not None:
            self._potential_energy = energy
            self.implemented_properties.append("energy")
        if charges is not None:
            self._charges = np.asfarray(charges)
            self.implemented_properties.append("charges")

    def _calculate(self, atoms: Atoms, **kwargs: Any) -> None:  # noqa: D102
        if "forces" in self.implemented_properties:
            self.results["forces"] = self._forces
        if "charges" in self.implemented_properties:
            self.results["charges"] = self._charges
        if "energy" in self.implemented_properties:
            self.results["energy"] = self._potential_energy

    def assign_atoms(self, atoms: Optional[Atoms]) -> None:  # noqa: D102
        if atoms is not None:
            if "forces" in self.implemented_properties and len(self._forces) != len(
                atoms
            ):
                raise CalculatorSetupError(
                    "`forces` array is not the same size as the system."
                )
            if "charges" in self.implemented_properties and len(self._charges) != len(
                atoms
            ):
                raise CalculatorSetupError(
                    "`charges` array is not the same size as the system."
                )
