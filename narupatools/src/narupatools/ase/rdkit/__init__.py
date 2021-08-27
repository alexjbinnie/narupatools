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

"""
Interoperability between ASE and RDKit.

This module includes code which is uses RDKit but through ASE. It exposes
both RDKit force fields (MMFF94(s) and UFF) as calculators that can be used
with ASE. Combined with the ability to quickly create ASE atoms objects from
SMILES strings, this allows quick generation of simulations.
"""

import importlib.util

if importlib.util.find_spec("rdkit") is None:
    raise ImportError("Cannot use narupatools.ase.rdkit without rdkit.")

from ._calculator import MMFF94Calculator, UFFCalculator
from ._generate import atoms_from_smiles

__all__ = ["atoms_from_smiles", "MMFF94Calculator", "UFFCalculator"]
