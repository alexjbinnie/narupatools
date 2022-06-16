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

pytest.importorskip("lammps")

from narupatools.lammps.exceptions import ComputeNotFoundError


def test_gather_compute_atom_1d(lammps):
    lammps.command("compute atom1 all property/atom x")
    value = lammps.gather_compute("atom1")
    assert value.shape == (2004,)
    assert value[0] == pytest.approx(43.99993)


def test_gather_compute_atom_2d(lammps):
    lammps.command("compute atom1 all property/atom x y z")
    value = lammps.gather_compute("atom1")
    assert value.shape == (2004, 3)
    assert value[0][0] == pytest.approx(43.99993)


def test_gather_atom_compute_missing(lammps):
    with pytest.raises(ComputeNotFoundError):
        lammps.gather_compute("atom1")
