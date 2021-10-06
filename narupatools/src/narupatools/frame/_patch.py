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

"""Patch to FrameData that allows copy() to work with empty arrays."""
from collections import KeysView
from typing import Any, Generator, Tuple, Union

import numpy as np
from narupa.trajectory import FrameData
from narupa.utilities.protobuf_utilities import value_to_object

from narupatools.frame.fields import FrameKey, get_frame_key


def _copy(self: FrameData) -> FrameData:
    frame = FrameData()
    for key, value in self.raw.arrays.items():
        frame.raw.arrays[key].CopyFrom(value)
    for key, value in self.raw.values.items():
        frame.raw.values[key].CopyFrom(value)
    return frame


def _repr(self: FrameData) -> str:
    rep = "<NarupaFrame"

    for key, value in self.items():
        rep += f" {key}={_print_value(value)}"

    rep += ">"

    return rep


def _keys(self) -> Generator[str, None, None]:
    """Iterate over the keys of the Frame."""
    yield from self.raw.values.keys()
    yield from self.raw.arrays.keys()


def _items(self) -> Generator[Tuple[str, Any], None, None]:
    """Iterate over the keys and values of the Frame."""
    for key in self.keys():
        yield key, self[key]


def _getitem(self, k: Union[str, FrameKey]) -> Any:
    if isinstance(k, FrameKey):
        return k.get(self)
    try:
        return get_frame_key(k).get(self)
    except KeyError as e:
        if k in self.raw.values:
            return value_to_object(self.raw.values[k])
        if k in self.raw.arrays:
            arr = self.raw.arrays[k]
            print(type(self.raw.arrays[k]))
            if self.raw.arrays[k].HasField("index_values"):
                return np.array(arr.ListFields()[0][1].values, dtype=int)
            elif self.raw.arrays[k].HasField("float_values"):
                return np.array(arr.ListFields()[0][1].values, dtype=float)
            elif self.raw.arrays[k].HasField("string_values"):
                return np.array(arr.ListFields()[0][1].values, dtype=object)
        raise KeyError from e

def _setitem(self, key: Union[str, FrameKey], value: Any) -> None:
    if isinstance(key, FrameKey):
        key.set(self, value)
        return
    try:
        get_frame_key(key).set(self, value)
    except KeyError:
        self._set_from_type(key, value)

def _set_from_type(self, key: str, value: Any) -> None:
    if isinstance(value, str):
        self.set_string_value(key, value)
    elif isinstance(value, float):
        self.set_float_value(key, value)
    elif isinstance(value, np.ndarray):
        if value.dtype == float:
            self.set_float_array(key, value)
        elif value.dtype == int:
            if all(i >= 0 for i in value):
                self.set_index_array(key, value)
            else:
                self.set_float_array(key, value)
        elif value.dtype == object:
            self.set_string_array(key, value)
        else:
            raise TypeError(f"Did not know how to serialize {value}.")
    else:
        raise TypeError(f"Did not know how to serialize {value}.")

def _contains(self, key: Any) -> bool:
    if isinstance(key, FrameKey):
        return key.key in self.arrays or key.key in self.values
    return key in self.arrays or key in self.values


def _print_value(value: Any) -> str:
    if isinstance(value, np.ndarray) and len(value) > 3:
        with np.printoptions(precision=3, suppress=True):
            return f"[{value[0]}, {value[1]}, ... ({len(value)} item(s)]"
    return repr(value)


FrameData.copy = _copy  # type: ignore
FrameData.__repr__ = _repr
FrameData.keys = _keys
FrameData.items = _items
FrameData.__getitem__ = _getitem
FrameData.__setitem__ = _setitem
FrameData.__contains__ = _contains
