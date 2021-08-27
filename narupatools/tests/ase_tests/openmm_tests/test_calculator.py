import copy

import pytest

from narupatools.ase.openmm import OpenMMCalculator


@pytest.fixture
def calculator(villin_openmm_simulation):
    return OpenMMCalculator(villin_openmm_simulation)


def test_calculate_specified_atoms(calculator, villin_ase_atoms):
    calculator.calculate(villin_ase_atoms)
