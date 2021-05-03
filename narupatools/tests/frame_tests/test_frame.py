from collections import Mapping

import numpy as np
import pytest
from narupa.trajectory import FrameData

from narupatools.frame import NarupaFrame


def test_frame():
    frame = NarupaFrame()
    assert isinstance(frame, Mapping)


def test_set_string_value():
    frame = NarupaFrame()
    frame["abc"] = "my_string"
    assert frame["abc"] == "my_string"


def test_set_float_value():
    frame = NarupaFrame()
    frame["abc"] = 2.4
    assert frame["abc"] == pytest.approx(2.4)


def test_set_float_array():
    frame = NarupaFrame()
    frame["abc"] = np.array([0.0, 1.0, 3.0], dtype=float)
    assert frame["abc"] == pytest.approx(np.array([0.0, 1.0, 3.0]))
    assert frame["abc"].dtype == float


def test_set_int_array_positive():
    frame = NarupaFrame()
    frame["abc"] = np.array([0, 1, 3], dtype=int)
    assert frame["abc"] == pytest.approx(np.array([0.0, 1.0, 3.0]))
    assert frame["abc"].dtype == int


def test_set_int_array_mixed():
    frame = NarupaFrame()
    frame["abc"] = np.array([-4, 1, 3], dtype=int)
    assert frame["abc"] == pytest.approx(np.array([-4.0, 1.0, 3.0]))
    assert frame["abc"].dtype == float


def test_set_string_array():
    frame = NarupaFrame()
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
    frame2 = NarupaFrame()
    frame2.raw.MergeFrom(frame.raw)
    assert frame2["array"] == [1.0]
