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
from ase.io import read
from testing import add_mark


def pytest_collection_modifyitems(items):
    add_mark(filename=__file__, mark=pytest.mark.ase, items=items)


@pytest.fixture(scope="session")
def ethane_atoms_readonly(ethane_sdf_filename) -> Atoms:
    return read(ethane_sdf_filename)  # type: ignore[return-value]


@pytest.fixture
def ethane_atoms(ethane_atoms_readonly) -> Atoms:
    return ethane_atoms_readonly.copy()  # type: ignore[no-any-return]


@pytest.fixture
def single_carbon_atoms() -> Atoms:
    return Atoms(
        symbols=["C"], masses=np.array([12.0]), positions=np.array([[0.0, 0.0, 0.0]])
    )
