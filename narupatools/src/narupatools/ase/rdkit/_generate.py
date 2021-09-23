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

from ase import Atoms
from rdkit import Chem
from rdkit.Chem import AllChem

from narupatools.frame import convert


def atoms_from_smiles(*smiles: str, add_hydrogens: bool = True) -> Atoms:
    """
    Generate an ASE Atoms object from a SMILES string.

    This leverages RDKit's ability to generate molecules from SMILES
    strings to create a system.

    :param smiles: SMILES string describing the system.
    :param add_hydrogens: Should implicit hydrogens be added to the system.
    :return: Atoms object representing system.
    """
    if not isinstance(smiles, str):
        smiles = '.'.join(smiles)
    mol = Chem.MolFromSmiles(smiles)
    if add_hydrogens:
        mol = AllChem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    # Minimize to avoid clashes
    AllChem.MMFFOptimizeMolecule(mol, ignoreInterfragInteractions=False)
    return convert(mol, Atoms)
