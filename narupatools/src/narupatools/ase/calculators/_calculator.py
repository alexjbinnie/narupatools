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

from __future__ import annotations

from abc import ABCMeta, abstractmethod
from collections import Collection
from typing import List, Optional, final

from ase import Atoms
from ase.calculators.calculator import Calculator as ASECalculator
from ase.calculators.calculator import CalculatorSetupError, all_changes

from narupatools.override import override


class Calculator(ASECalculator, metaclass=ABCMeta):
    """
    Subclass of ASE calculator which has additional callbacks.

    This calculator has two features. First, it has a method :meth:`assign_atoms` which is called
    when the atoms object (which represents the last atoms object used in a calculation) is set.
    This function can be used by subclasses to check the atoms object is valid.

    This calculator also implements the :meth:`calculate` function that all calculators require. It
    performs some general logic such as checking that there is actually an atoms object to use for
    a calculation. Subclasses should instead override :meth:`_calculate`, which instead of an
    Optional[Atoms] object takes just an Atoms object. This means subclasses do not have to worry
    about checks if the atoms object is None when it comes to typing, and they should use this
    argument instead of self.atoms (which may be None).
    """

    @property  # type: ignore[override]
    def atoms(self) -> Optional[Atoms]:  # type: ignore[override] # noqa: D102
        return self._atoms

    @atoms.setter
    def atoms(self, atoms: Optional[Atoms]) -> None:
        self._atoms = atoms
        self.assign_atoms(atoms)

    def assign_atoms(self, atoms: Optional[Atoms]) -> None:
        """
        Perform additional checks when an :class:`Atoms` object is assigned to this calculator.

        This assignment may occur at a couple of points, including initial initialization of the
        calculator, or before :meth:`calculate` is called if the Atoms object has changed. As for
        the base class, this Atoms object is a copy of the one being used in the calculator.

        :param atoms: Atoms object which has been assigned to this calculator.
        """

    @override(ASECalculator.calculate)
    @final
    def calculate(  # noqa: D102
        self,
        atoms: Optional[Atoms] = None,
        properties: Collection[str] = ("forces", "energy"),
        system_changes: List[str] = all_changes,
    ) -> None:
        super().calculate(atoms, properties, system_changes)
        if self.atoms is None:
            raise CalculatorSetupError("Atoms not provided to calculator.")
        self._calculate(
            self.atoms, properties=properties, system_changes=system_changes
        )

    @abstractmethod
    def _calculate(
        self,
        atoms: Atoms,
        *,
        properties: Collection[str] = ("forces", "energy"),
        system_changes: List[str] = all_changes
    ) -> None:
        pass
