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

import abc
from typing import Any, List, Optional

from ase import Atoms
from ase.calculators.calculator import all_changes
from rdkit import Chem

from narupatools.ase import UnitsASE
from narupatools.core import UnitsNarupa
from narupatools.frame import convert

from ..calculators import Calculator
from ._forcefield import _MMFFForceField, _RDKitForcefield, _UFFForceField

_ASEToNarupa = UnitsASE >> UnitsNarupa
_NarupaToASE = UnitsNarupa >> UnitsASE

_ASE_TO_NARUPA_LENGTH = _ASEToNarupa.length
_NARUPA_TO_ASE_FORCE = _NarupaToASE.force
_NARUPA_TO_ASE_ENERGY = _NarupaToASE.energy


class RDKitForceFieldCalculator(Calculator, metaclass=abc.ABCMeta):
    """
    ASE calculator that performs calculations using a forcefield included within RDKit.

    There are two calculators, :class:`MMFF94Calculator` and :class:`UFFCalculator`, which
    can be used. These calculators will throw an error if parameters are not found for all
    the atoms in the system.

    Note that both these forcefields require bonds. Bond information is added to ASE by
    narupatools in some cases.
    """

    implemented_properties = ["energy", "forces"]

    def __init__(self, **kwargs: Any):
        super().__init__(**kwargs)
        self._forcefield: Optional[_RDKitForcefield] = None

    @abc.abstractmethod
    def _create_forcefield(self, atoms: Atoms) -> _RDKitForcefield:
        ...

    def _calculate(
        self, atoms: Atoms, *, system_changes: List[str] = all_changes, **kwargs: Any
    ) -> None:
        if "numbers" in system_changes or self._forcefield is None:
            self._forcefield = self._create_forcefield(atoms)

        if "numbers" in system_changes or "positions" in system_changes:
            self._calculate_forcefield(atoms)

    def _calculate_forcefield(self, atoms: Atoms) -> None:

        if self._forcefield is None:
            raise ValueError("Forcefield is null.")
        self.results["forces"] = (
            self._forcefield.calculate_forces(atoms.positions * _ASE_TO_NARUPA_LENGTH)
            * _NARUPA_TO_ASE_FORCE
        )
        self.results["energy"] = (
            self._forcefield.calculate_energy(atoms.positions * _ASE_TO_NARUPA_LENGTH)
            * _NARUPA_TO_ASE_ENERGY
        )


class UFFCalculator(RDKitForceFieldCalculator):
    """ASE calculator that used the UFF force field in RDKit."""

    default_parameters = {"nonbonded_cutoff": 1.0, "include_interatomic": True}

    def _create_forcefield(self, atoms: Atoms) -> _UFFForceField:
        mol = convert(atoms, Chem.rdchem.Mol)
        return _UFFForceField(
            mol,
            nonbonded_cutoff=self.parameters.nonbonded_cutoff,
            include_interatomic=self.parameters.include_interatomic,
        )


class MMFF94Calculator(RDKitForceFieldCalculator):
    """ASE calculator that used the UFF force field in RDKit."""

    default_parameters = {
        "variant": "mmff94",
        "nonbonded_cutoff": 10.0,
        "include_interatomic": True,
    }

    def _create_forcefield(self, atoms: Atoms) -> _MMFFForceField:
        mol = convert(atoms, Chem.rdchem.Mol)
        return _MMFFForceField(
            mol,
            variant=self.parameters.variant,
            nonbonded_cutoff=self.parameters.nonbonded_cutoff,
            include_interatomic=self.parameters.include_interatomic,
        )
