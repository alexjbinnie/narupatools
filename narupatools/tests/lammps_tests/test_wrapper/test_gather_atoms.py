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

from narupatools.lammps.exceptions import UnknownAtomPropertyError


def test_gather_atoms_float_3d(lammps):
    value = lammps.gather_atoms("f", dimension=3)
    assert isinstance(value, np.ndarray)
    assert value.shape == (2004, 3)


def test_gather_atoms_float_1d(lammps):
    value = lammps.gather_atoms("q", dimension=1)
    assert isinstance(value, np.ndarray)
    assert value.shape == (2004,)


def test_gather_atoms_int_1d(lammps):
    value = lammps.gather_atoms("id", dimension=1)
    assert isinstance(value, np.ndarray)
    assert value.shape == (2004,)
    assert value.dtype == np.int32


def test_gather_atoms_invalid(lammps):
    with pytest.raises(UnknownAtomPropertyError):
        lammps.gather_atoms("blah", 3)


def test_gather_atoms_compute(lammps):
    lammps.command("compute atom1 all property/atom x")
    with pytest.raises(UnknownAtomPropertyError):
        lammps.gather_atoms("c_atom1", dimension=1)


def test_positions(lammps):
    value = lammps.gather_atoms("x", dimension=3)
    assert value.shape == (2004, 3)
    assert value[0] == pytest.approx(np.array([43.99993, 58.52678, 36.7855]))
