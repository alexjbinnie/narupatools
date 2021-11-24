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

"""Code for interface with the RDKit package."""

import importlib

has_rdkit = importlib.util.find_spec("rdkit") is not None

if not has_rdkit:
    raise ImportError("narupatools.rdkit requires rdkit to be installed.")

from ._converter import frame_to_rdkit_mol, rdkit_mol_to_frame
from ._generate import generate_from_smiles
from ._units import UnitsRDKit

__all__ = [
    "generate_from_smiles",
    "frame_to_rdkit_mol",
    "rdkit_mol_to_frame",
    "UnitsRDKit",
]
