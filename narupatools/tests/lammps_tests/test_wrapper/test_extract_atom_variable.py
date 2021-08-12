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

from narupatools.lammps.exceptions import VariableNotFoundError


def test_extract_atom_variable(lammps):
    lammps.command("variable my_var atom x*2.0")
    value = lammps.extract_atom_variable("my_var")
    assert isinstance(value, np.ndarray)
    assert value.shape == (2004,)
    assert value[0] == pytest.approx(78.01242455, rel=1e-3)


def test_extract_atom_variable_missing(lammps):
    with pytest.raises(VariableNotFoundError):
        lammps.extract_atom_variable("my_var")


def test_extract_atom_variable_none_key(lammps):
    with pytest.raises(VariableNotFoundError):
        lammps.extract_atom_variable(None)


def test_extract_atom_variable_int_key(lammps):
    with pytest.raises(VariableNotFoundError):
        lammps.extract_atom_variable(3.1)
