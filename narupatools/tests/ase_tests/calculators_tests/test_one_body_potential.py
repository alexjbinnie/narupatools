import numpy as np
import pytest
from ase import Atom

from narupatools.ase import OneBodyPotentialCalculator
from narupatools.physics.typing import Vector3
from narupatools.physics.vector import magnitude


class HarmonicWellPotentialCalculator(OneBodyPotentialCalculator):
    def calculate_energy(self, atom: Atom) -> float:
        return magnitude(atom.position) ** 2

    def calculate_force(self, atom: Atom) -> Vector3:
        return -2.0 * np.array(atom.position)


@pytest.fixture
def carbon_atoms(single_carbon_atoms):
    single_carbon_atoms.positions[0] = [2.0, 0.0, 0.0]
    return single_carbon_atoms


@pytest.fixture
def well_calculator():
    return HarmonicWellPotentialCalculator()


def test_has_forces(carbon_atoms, well_calculator):
    well_calculator.calculate(carbon_atoms)
    assert well_calculator.results["forces"] == pytest.approx(
        np.array([[-4.0, 0.0, 0.0]])
    )


def test_has_energies(carbon_atoms, well_calculator):
    well_calculator.calculate(carbon_atoms)
    assert well_calculator.results["energy"] == pytest.approx(4.0)


def test_get_forces(carbon_atoms, well_calculator):
    carbon_atoms.set_calculator(well_calculator)
    assert carbon_atoms.get_forces() == pytest.approx(np.array([[-4.0, 0.0, 0.0]]))


def test_get_potential_energy(carbon_atoms, well_calculator):
    carbon_atoms.set_calculator(well_calculator)
    assert carbon_atoms.get_potential_energy() == pytest.approx(4.0)
