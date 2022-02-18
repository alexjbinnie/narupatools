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

from narupatools.physics.matrix import kronecker_delta, zero_matrix


def test_zero_matrix():
    assert zero_matrix() == pytest.approx(
        np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    )


@pytest.mark.parametrize("i", [0, 1, 3])
@pytest.mark.parametrize("j", [0, 1, 3])
def test_kronecker_delta(i, j):
    assert kronecker_delta(i, j) == (1 if i == j else 0)
