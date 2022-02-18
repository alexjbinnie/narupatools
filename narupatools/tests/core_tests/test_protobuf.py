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
from narupa.utilities.protobuf_utilities import dict_to_struct

from narupatools.core.protobuf import dictionary_to_protobuf


def test_set():
    dict_ = {"my_set": {1, 3, 5}}
    _ = dictionary_to_protobuf(dict_)
    with pytest.raises(ValueError):  # noqa: PT011
        _ = dict_to_struct(dict_)


def test_range():
    dict_ = {"my_set": range(10)}
    _ = dictionary_to_protobuf(dict_)
    with pytest.raises(ValueError):  # noqa: PT011
        _ = dict_to_struct(dict_)


def test_numpy():
    dict_ = {"my_set": np.array([0.0, 1.0])}
    _ = dictionary_to_protobuf(dict_)
    with pytest.raises(ValueError):  # noqa: PT011
        _ = dict_to_struct(dict_)
