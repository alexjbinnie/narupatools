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

"""Protocols for use with ASE calculators."""

from typing import Protocol

from ase import Atoms


class CalculatorSetAtoms(Protocol):
    """
    Protocol for an ASE calculator that defines the set_atoms method.

    This protocol is a useful base class to ensure a calculator implements this method
    correctly, and provides a docstring describing it.
    """

    def set_atoms(self, atoms: Atoms, /) -> None:
        """
        Called when an :class:`ase.atoms.Atoms` object is assigned this calculator via atoms.calc.

        :param atoms: Atoms this is assigned to.
        """
        ...
