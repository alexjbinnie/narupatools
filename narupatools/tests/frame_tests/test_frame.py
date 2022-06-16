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
from narupa.trajectory import FrameData


def test_set_string_value():
    frame = FrameData()
    frame["abc"] = "my_string"
    assert frame["abc"] == "my_string"


def test_set_float_value():
    frame = FrameData()
    frame["abc"] = 2.4
    assert frame["abc"] == pytest.approx(2.4)


def test_set_float_array():
    frame = FrameData()
    frame["abc"] = np.array([0.0, 1.0, 3.0], dtype=float)
    assert frame["abc"] == pytest.approx(np.array([0.0, 1.0, 3.0]))
    assert frame["abc"].dtype == float


def test_set_int_array_positive():
    frame = FrameData()
    frame["abc"] = np.array([0, 1, 3], dtype=int)
    assert frame["abc"] == pytest.approx(np.array([0.0, 1.0, 3.0]))
    assert frame["abc"].dtype == int


def test_set_int_array_mixed():
    frame = FrameData()
    frame["abc"] = np.array([-4, 1, 3], dtype=int)
    assert frame["abc"] == pytest.approx(np.array([-4.0, 1.0, 3.0]))
    assert frame["abc"].dtype == float


def test_set_string_array():
    frame = FrameData()
    frame["abc"] = np.array(["a", "b", "c"], dtype=object)
    assert np.all(frame["abc"] == np.array(["a", "b", "c"]))
    assert frame["abc"].dtype == object


def test_copy_empty_array():
    frame = FrameData()
    frame.set_float_array("array", [])
    frame2 = frame.copy()
    assert "array" in frame2.arrays


def test_copy_non_empty():
    frame = FrameData()
    frame.set_float_array("array", [1.0, 0.0])
    frame2 = frame.copy()
    assert frame2.arrays["array"] == pytest.approx([1.0, 0.0])


def test_merge_from():
    frame = FrameData()
    frame.set_float_array("array", [1.0])
    frame2 = FrameData()
    frame2.raw.MergeFrom(frame.raw)
    assert frame2["array"] == [1.0]
