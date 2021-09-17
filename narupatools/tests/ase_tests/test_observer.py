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
from testing import assert_event_called

from narupatools.ase.constraints._observer import ASEObserver


@pytest.fixture
def observer(ethane_atoms):
    return ASEObserver.get_or_create(ethane_atoms)


def test_set_positions(ethane_atoms, observer):
    positions = ethane_atoms.get_positions()
    with assert_event_called(observer.on_set_positions):
        ethane_atoms.set_positions(positions)


def test_set_momenta(ethane_atoms, observer):
    momenta = ethane_atoms.get_momenta()
    with assert_event_called(observer.on_set_momenta):
        ethane_atoms.set_momenta(momenta)


def test_set_cell(ethane_atoms, observer):
    cell = ethane_atoms.get_cell()
    with assert_event_called(observer.on_set_cell):
        ethane_atoms.set_cell(cell)
