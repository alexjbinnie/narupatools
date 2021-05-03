import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import CalculatorSetupError, PropertyNotImplementedError
from ase.data.g2 import data

from narupatools.ase.calculators import ConstantCalculator
from narupatools.physics.vector import vector


@pytest.fixture
def benzene():
    mol = data["C6H6"]
    return Atoms(positions=mol["positions"], symbols=mol["symbols"])


def test_atoms_arg(benzene):
    calc = ConstantCalculator(atoms=benzene)
    assert benzene.calc is calc


def test_no_args(benzene):
    calc = ConstantCalculator()
    benzene.calc = calc
    with pytest.raises(PropertyNotImplementedError):
        benzene.get_forces()
    with pytest.raises(PropertyNotImplementedError):
        benzene.get_potential_energy()
    with pytest.raises(PropertyNotImplementedError):
        benzene.get_charges()


def test_valid(benzene):
    forces = [
        vector(1, 0, 0),
        vector(0, 1, 0),
        vector(1, 0, 0),
        vector(2, 0, 0),
        vector(0, 1, 0),
        vector(-1, 0, 0),
        vector(0, 3, -1),
        vector(1, 1, 0),
        vector(0, 0, 1),
        vector(0, 0, 2),
        vector(3, 0, -1),
        vector(0, 2, 1),
    ]
    charges = [0, 2, 5, -2, 1, 0, 0, 3, -2, 0, 2, 0]
    energy = -23.21
    calc = ConstantCalculator(forces=forces, charges=charges, energy=energy)
    benzene.calc = calc
    assert benzene.get_forces() == pytest.approx(np.asfarray(forces))
    assert benzene.get_potential_energy() == pytest.approx(energy)
    assert benzene.get_charges() == pytest.approx(charges)


def test_forces_invalid(benzene):
    forces = [
        vector(1, 0, 0),
        vector(0, 1, 0),
        vector(1, 0, 0),
        vector(2, 0, 0),
        vector(0, 1, 0),
        vector(-1, 0, 0),
        vector(0, 3, -1),
        vector(1, 1, 0),
        vector(0, 0, 1),
        vector(0, 0, 2),
    ]
    calc = ConstantCalculator(forces=forces)
    with pytest.raises(CalculatorSetupError):
        benzene.calc = calc


def test_charges_invalid(benzene):
    charges = [0, 2, 5, -2, 1, 0, 0, 2, 0]
    calc = ConstantCalculator(charges=charges)
    with pytest.raises(CalculatorSetupError):
        benzene.calc = calc
