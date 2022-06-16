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

from narupatools.lammps.exceptions import VariableNotFoundError


def test_extract_variable(lammps):
    lammps.command("variable my_var equal temp/3.0")
    value = lammps.extract_variable("my_var")
    assert isinstance(value, float)
    assert value == pytest.approx(94.03350571744606, rel=1e-3)


def test_extract_variable_missing(lammps):
    with pytest.raises(VariableNotFoundError):
        lammps.extract_variable("my_var")
