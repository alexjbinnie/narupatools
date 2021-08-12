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

from narupatools.lammps.exceptions import ComputeNotFoundError


def test_extract_local_compute_1d(lammps):
    lammps.command("compute bond1 all property/local btype")
    value = lammps.extract_local_compute("bond1")
    assert value.shape == (1365,)
    assert not value.flags.writeable


def test_extract_local_compute_2d(lammps):
    lammps.command("compute bond1 all property/local btype batom1 batom2")
    value = lammps.extract_local_compute("bond1")
    assert value.shape == (1365, 3)
    assert not value.flags.writeable


def test_extract_local_compute_missing(lammps):
    with pytest.raises(ComputeNotFoundError):
        lammps.extract_local_compute("bond1")
