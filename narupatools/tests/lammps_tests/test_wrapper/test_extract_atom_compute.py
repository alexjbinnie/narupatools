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


def test_extract_atom_compute_1d(lammps):
    lammps.command("compute atom1 all property/atom x")
    value = lammps.extract_atom_compute("atom1")
    assert value.shape == (2004,)
    assert not value.flags.writeable


def test_extract_atom_compute_2d(lammps):
    lammps.command("compute atom1 all property/atom x y z fx fy fz")
    value = lammps.extract_atom_compute("atom1")
    assert value.shape == (2004, 6)
    assert not value.flags.writeable


def test_extract_atom_compute_ke(lammps):
    lammps.command("compute atom_ke all ke/atom")
    value = lammps.extract_atom_compute("atom_ke")
    assert value.shape == (2004,)
    assert not value.flags.writeable


def test_extract_atom_compute_missing(lammps):
    with pytest.raises(ComputeNotFoundError):
        lammps.extract_atom_compute("atom1")
