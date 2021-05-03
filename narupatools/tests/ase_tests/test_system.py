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

import pytest
from infinite_sets import everything

from narupatools.ase import NullCalculator
from narupatools.ase.system import ASESystem
from narupatools.frame import BondPairs


@pytest.fixture(scope="module")
def system(villin_ase_atoms_readonly) -> ASESystem:
    atoms = villin_ase_atoms_readonly.copy()
    atoms.set_calculator(NullCalculator())
    return ASESystem(atoms)


@pytest.fixture(scope="module")
def system_universe(
    villin_ase_atoms_readonly, villin_mda_universe_readonly
) -> ASESystem:
    atoms = villin_ase_atoms_readonly.copy()
    atoms.set_calculator(NullCalculator())
    return ASESystem(atoms, universe=villin_mda_universe_readonly)


def test_system_missing_bonds(system):
    frame = system.get_frame(everything())
    assert BondPairs.key not in frame


def test_system_universe_has_bonds(system_universe):
    frame = system_universe.get_frame(everything())
    assert BondPairs.key in frame
