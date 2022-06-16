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
