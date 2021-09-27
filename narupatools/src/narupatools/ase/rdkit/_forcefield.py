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

from typing import Any, Literal, Optional

import numpy as np
from ase.calculators.calculator import CalculatorSetupError
from rdkit import Chem
from rdkit.Chem.rdForceFieldHelpers import (
    MMFFGetMoleculeForceField,
    MMFFGetMoleculeProperties,
    MMFFHasAllMoleculeParams,
    UFFGetMoleculeForceField,
    UFFHasAllMoleculeParams,
)

from narupatools.physics.units import UnitsNarupa
from narupatools.rdkit import UnitsRDKit

_RDKitToNarupa = UnitsRDKit >> UnitsNarupa
_NarupaToRDKit = UnitsNarupa >> UnitsRDKit


class _RDKitForcefield:
    def __init__(self, forcefield: Any):
        self._forcefield = forcefield

    def calculate_energy(self, coordinates: Optional[np.ndarray] = None) -> float:
        if coordinates is not None:
            energy = self._forcefield.CalcEnergy(
                list(coordinates.flatten() * _NarupaToRDKit.length)
            )
        else:
            energy = self._forcefield.CalcEnergy()
        return energy * _RDKitToNarupa.energy  # type: ignore[no-any-return]

    def calculate_forces(self, coordinates: Optional[np.ndarray] = None) -> np.ndarray:
        if coordinates is not None:
            forces = self._forcefield.CalcGrad(
                list(coordinates.flatten() * _NarupaToRDKit.length)
            )
        else:
            forces = self._forcefield.CalcGrad()
        return -np.array(forces).reshape((-1, 3)) * _RDKitToNarupa.force  # type: ignore[no-any-return]


class _MMFFForceField(_RDKitForcefield):
    """Narupa-friendly wrapper around the MMFF94(s) forcefield as implemented in RDKit."""

    def __init__(
        self,
        mol: Chem.rdchem.Mol,
        nonbonded_cutoff: float = 10,
        include_interatomic: bool = True,
        variant: Literal["mmff94", "mmff94s"] = "mmff94",
    ):
        """
        Create an MMFF force field.

        :param mol: RDKit molecule to parameterize.
        :param nonbonded_cutoff: Nonbonded cutoff in nanometers.
        :param include_interatomic: Should forces between separate fragments be included.
        :param variant: What variant of mmff94 to use.
        """
        has_params = MMFFHasAllMoleculeParams(mol)
        if not has_params:
            raise CalculatorSetupError("Cannot find MMFF parameters.")
        properties = MMFFGetMoleculeProperties(mol, mmffVariant=variant)
        force_field = MMFFGetMoleculeForceField(
            mol,
            properties,
            nonBondedThresh=nonbonded_cutoff,
            ignoreInterfragInteractions=not include_interatomic,
        )
        super().__init__(force_field)


class _UFFForceField(_RDKitForcefield):
    """Narupa-friendly wrapper around the UFF forcefield as implemented in RDKit."""

    def __init__(
        self,
        mol: Chem.rdchem.Mol,
        nonbonded_cutoff: float = 1,
        include_interatomic: bool = True,
    ):
        """
        Create an UFF force field.

        :param mol: RDKit molecule to parameterize.
        :param nonbonded_cutoff: Nonbonded cutoff in nanometers.
        :param include_interatomic: Should forces between separate fragments be included.
        """
        has_params = UFFHasAllMoleculeParams(mol)
        if not has_params:
            raise CalculatorSetupError("Cannot find UFF parameters.")
        force_field = UFFGetMoleculeForceField(
            mol,
            nonbonded_cutoff * _NarupaToRDKit.length,
            ignoreInterfragInteractions=not include_interatomic,
        )
        super().__init__(force_field)
