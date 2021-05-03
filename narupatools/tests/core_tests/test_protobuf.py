import numpy as np
import pytest
from narupa.utilities.protobuf_utilities import dict_to_struct

from narupatools.core.protobuf import dictionary_to_protobuf


def test_set():
    dict = {"my_set": {1, 3, 5}}
    _ = dictionary_to_protobuf(dict)
    with pytest.raises(ValueError):  # noqa: PT011
        _ = dict_to_struct(dict)


def test_range():
    dict = {"my_set": range(10)}
    _ = dictionary_to_protobuf(dict)
    with pytest.raises(ValueError):  # noqa: PT011
        _ = dict_to_struct(dict)


def test_numpy():
    dict = {"my_set": np.array([0.0, 1.0])}
    _ = dictionary_to_protobuf(dict)
    with pytest.raises(ValueError):  # noqa: PT011
        _ = dict_to_struct(dict)
