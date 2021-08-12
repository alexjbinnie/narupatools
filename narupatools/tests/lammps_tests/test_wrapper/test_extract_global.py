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

from narupatools.lammps.exceptions import GlobalNotFoundError


def test_extract_global_int(lammps):
    value = lammps.extract_global("nlocal")
    assert isinstance(value, int)
    assert value == 2004


def test_extract_global_float(lammps):
    value = lammps.extract_global("dt")
    assert isinstance(value, float)
    assert value == pytest.approx(2.0)


def test_extract_global_string(lammps):
    value = lammps.extract_global("units")
    assert isinstance(value, str)
    assert value == "real"


def test_extract_global_int_autodetect(lammps):
    value = lammps.extract_global("nlocal")
    assert isinstance(value, int)
    assert value == 2004


def test_extract_global_float_autodetect(lammps):
    value = lammps.extract_global("dt")
    assert isinstance(value, float)
    assert value == pytest.approx(2.0)


def test_extract_global_string_autodetect(lammps):
    value = lammps.extract_global("units")
    assert isinstance(value, str)
    assert value == "real"


def test_extract_global_missing_autodetect(lammps):
    with pytest.raises(GlobalNotFoundError):
        lammps.extract_global("missing_global")


def test_extract_global_missing_int(lammps):
    with pytest.raises(GlobalNotFoundError):
        lammps.extract_global("missing_global")


def test_extract_global_missing_float(lammps):
    with pytest.raises(GlobalNotFoundError):
        lammps.extract_global("missing_global")


def test_extract_global_missing_string(lammps):
    with pytest.raises(GlobalNotFoundError):
        lammps.extract_global("missing_global")
