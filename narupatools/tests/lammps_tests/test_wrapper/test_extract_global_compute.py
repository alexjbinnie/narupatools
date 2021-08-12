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

from narupatools.lammps.constants import VariableDimension
from narupatools.lammps.exceptions import ComputeNotFoundError


def test_built_in_compute(lammps):
    value = lammps.extract_global_compute(
        "thermo_pe", dimension=VariableDimension.SCALAR
    )
    assert isinstance(value, float)


def test_compute_vector(lammps):
    lammps.command("compute my_compute all temp/com")
    vector = lammps.extract_global_compute(
        "my_compute", dimension=VariableDimension.VECTOR1D
    )
    assert isinstance(vector, np.ndarray)
    assert vector.shape == (6,)
    scalar = lammps.extract_global_compute(
        "my_compute", dimension=VariableDimension.SCALAR
    )
    assert isinstance(scalar, float)


def test_extract_compute_missing(lammps):
    with pytest.raises(ComputeNotFoundError):
        lammps.extract_global_compute("missing_key", dimension=VariableDimension.SCALAR)
