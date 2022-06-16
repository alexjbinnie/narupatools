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

from typing import Optional

import pytest
from ase import Atoms
from ase.calculators.calculator import all_changes

from narupatools.ase.calculators._calculator import Calculator


class CalculateOnPositionChangedCalculator(Calculator):

    calculation_count = 0
    assign_count = 0
    implemented_properties = ["energy"]

    def _calculate(self, atoms, system_changes=all_changes, **kwargs):
        if "positions" in system_changes:
            self.calculation_count += 1
        self.results["energy"] = 0.0

    def assign_atoms(self, atoms: Optional[Atoms]):
        if atoms is not None:
            self.assign_count += 1


class CalculateForeverCalculator(Calculator):

    calculation_count = 0
    assign_count = 0
    implemented_properties = ["energy"]
    ignored_changes = all_changes

    def _calculate(self, atoms, system_changes=all_changes, **kwargs):
        if "positions" in system_changes:
            self.calculation_count += 1
        self.results["energy"] = 0.0

    def assign_atoms(self, atoms: Optional[Atoms]):
        if atoms is not None:
            self.assign_count += 1


@pytest.fixture
def atom():
    return Atoms(symbols=["H"], masses=[1.0], positions=[[0.0, 0.0, 0.0]])


def test_calculator_position_changed(atom):
    atom.calc = CalculateOnPositionChangedCalculator()
    assert atom.calc.assign_count == 0
    assert atom.calc.calculation_count == 0
    _ = atom.get_potential_energy()
    assert atom.calc.assign_count == 1
    assert atom.calc.calculation_count == 1
    _ = atom.get_potential_energy()
    assert atom.calc.assign_count == 1
    assert atom.calc.calculation_count == 1
    atom.positions += [[1.0, 0.0, 0.0]]
    _ = atom.get_potential_energy()
    assert atom.calc.assign_count == 2
    assert atom.calc.calculation_count == 2


def test_calculator_ignore_changed(atom):
    atom.calc = CalculateForeverCalculator()
    assert atom.calc.assign_count == 0
    assert atom.calc.calculation_count == 0
    _ = atom.get_potential_energy()
    assert atom.calc.assign_count == 1
    assert atom.calc.calculation_count == 1
    _ = atom.get_potential_energy()
    assert atom.calc.assign_count == 1
    assert atom.calc.calculation_count == 1
    atom.positions += [[1.0, 0.0, 0.0]]
    _ = atom.get_potential_energy()
    assert atom.calc.assign_count == 1
    assert atom.calc.calculation_count == 1
    atom.calc.reset()
    _ = atom.get_potential_energy()
    assert atom.calc.assign_count == 2
    assert atom.calc.calculation_count == 2
