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

from narupatools.ase import NullCalculator
from narupatools.ase import UnitsASE
from narupatools.ase.constraints import InteractionConstraint
from narupatools.core import UnitsNarupa
from narupatools.imd import constant_interaction

ASEToNarupa = UnitsASE >> UnitsNarupa


@pytest.fixture
def carbon_atoms(single_carbon_atoms):
    single_carbon_atoms.positions[0] = [2.0, 0.0, 0.0]
    single_carbon_atoms.set_calculator(NullCalculator())
    return single_carbon_atoms


@pytest.fixture
def constraint(carbon_atoms):
    interaction = constant_interaction(force=[-1.0, 0.0, 0.0], particles=[0])
    constraint = InteractionConstraint(
        dynamics=carbon_atoms,
        key="abc",
        interaction=interaction,
        start_time=0.0,
    )
    carbon_atoms.constraints.append(constraint)
    return constraint


def test_constraint_force(carbon_atoms, constraint):
    assert carbon_atoms.get_forces() * ASEToNarupa.force == pytest.approx(
        np.array([[-1.0, 0.0, 0.0]])
    )


def test_constraint_energy(carbon_atoms, constraint):
    assert carbon_atoms.get_potential_energy() * ASEToNarupa.energy == pytest.approx(
        0.2
    )
