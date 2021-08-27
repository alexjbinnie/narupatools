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
