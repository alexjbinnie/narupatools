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
from ase import Atoms
from ase.calculators.calculator import CalculatorSetupError, PropertyNotImplementedError
from ase.calculators.mixing import AverageCalculator
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


def test_energy_avg(benzene):
    calc1 = ConstantCalculator(energy=10.0)
    calc2 = ConstantCalculator(energy=20.0)
    AverageCalculator([calc1, calc2], benzene)

    assert benzene.get_potential_energy() == pytest.approx(15.0)


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
    benzene.calc = ConstantCalculator(forces=forces)
    with pytest.raises(CalculatorSetupError):
        benzene.get_forces()


def test_charges_invalid(benzene):
    charges = [0, 2, 5, -2, 1, 0, 0, 2, 0]
    benzene.calc = ConstantCalculator(charges=charges)
    with pytest.raises(CalculatorSetupError):
        benzene.get_charges()
