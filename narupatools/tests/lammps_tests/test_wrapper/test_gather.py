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


def test_gather_fix_1d(lammps):
    lammps.command("fix my_fix all ave/atom 1 10 10 vx")
    lammps.command("run 100")
    value = lammps.gather("f_my_fix", dimension=1)
    assert isinstance(value, np.ndarray)
    assert value.shape == (2004,)
    assert value.dtype == np.float64
    assert not value.flags.writeable


def test_gather_fix_3d(lammps):
    lammps.command("fix my_fix all ave/atom 1 10 10 vx vy vz")
    lammps.command("run 100")
    value = lammps.gather("f_my_fix", dimension=3)
    assert isinstance(value, np.ndarray)
    assert value.shape == (2004, 3)
    assert value.dtype == np.float64
    assert not value.flags.writeable


def test_gather_compute_1d(lammps):
    lammps.command("compute atom1 all property/atom x")
    value = lammps.gather("c_atom1", dimension=1)
    assert value.shape == (2004,)
    assert value[0] == pytest.approx(43.99993)


def test_gather_compute_2d(lammps):
    lammps.command("compute atom1 all property/atom x y z")
    value = lammps.gather("c_atom1", dimension=3)
    assert value.shape == (2004, 3)
    assert value[0][0] == pytest.approx(43.99993)
