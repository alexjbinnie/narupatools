import pytest
from ase.calculators.calculator import CalculatorSetupError

from narupatools.ase.openmm import OpenMMCalculator


@pytest.fixture
def calculator(villin_openmm_simulation):
    return OpenMMCalculator(villin_openmm_simulation)


def test_calculate_no_atoms_raises_valueerror(calculator):
    with pytest.raises(CalculatorSetupError):
        calculator.calculate()


def test_calculate_specified_atoms(calculator, villin_ase_atoms):
    calculator.calculate(villin_ase_atoms)
