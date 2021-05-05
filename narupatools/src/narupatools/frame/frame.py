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

"""Extension of the existing Narupa FrameData."""

from __future__ import annotations

from collections.abc import Mapping
from typing import (
    Any,
    Generator,
    Iterator,
    Optional,
    SupportsFloat,
    Tuple,
    TypeVar,
    Union,
    overload,
)

import numpy as np
from narupa.trajectory import FrameData
from narupa.utilities.protobuf_utilities import value_to_object

from narupatools.frame.fields import FrameKey, get_frame_key

TFrom = TypeVar("TFrom")
TTo = TypeVar("TTo")
TDefault = TypeVar("TDefault")


class NarupaFrame(FrameData):
    """
    Wrapper around a protobuf FrameData.

    This subclass explicitly is an implementation of a mutable mapping between string
    keys and arbitrary values.
    """

    @overload
    def __setitem__(self, key: str, value: Any) -> None:
        ...

    @overload
    def __setitem__(self, key: FrameKey[TFrom, TTo], value: TFrom) -> None:
        ...

    def __setitem__(self, key: Union[str, FrameKey], value: Any) -> None:
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

    @overload
    def __getitem__(self, k: str) -> Any:
        ...

    @overload
    def __getitem__(self, k: FrameKey[TFrom, TTo]) -> TTo:
        ...

    def __getitem__(self, k: Union[str, FrameKey]) -> Any:
        if isinstance(k, FrameKey):
            return k.get(self)
        try:
            return get_frame_key(k).get(self)
        except KeyError:
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
            raise KeyError

    def __len__(self) -> int:
        return len(self.raw.values) + len(self.raw.arrays)

    def __iter__(self) -> Iterator[Any]:
        for key in self.keys():
            yield self[key]

    def set_string_value(self, key: str, value: Any) -> None:
        """
        Set a single string value with the given key.

        :param key: Key to set.
        :param value: Value to set.
        """
        self.raw.values[key].string_value = str(value)

    def set_float_value(self, key: str, value: SupportsFloat) -> None:
        """
        Set a single float value with the given key.

        :param key: Key to set.
        :param value: Value to set.
        """
        self.raw.values[key].number_value = float(value)

    def get_string_value(self, key: str) -> str:
        """
        Get a string value with a given key.

        :param key: Key to lookup.
        :return: String value for the given key.
        :raises KeyError: Given key is absent in this frame.
        """
        return self.raw.values[key].string_value  # type: ignore

    def get_float_value(self, key: str) -> float:
        """
        Get a float value with a given key.

        :param key: Key to lookup.
        :return: String value for the given key.
        :raises KeyError: Given key is absent in this frame.
        """
        return self.raw.values[key].number_value  # type: ignore

    def keys(self) -> Generator[str, None, None]:
        """Iterate over the keys of the Frame."""
        yield from self.raw.values.keys()
        yield from self.raw.arrays.keys()

    def items(self) -> Generator[Tuple[str, Any], None, None]:
        """Iterate over the keys and values of the Frame."""
        for key in self.keys():
            yield key, self[key]

    def values(self) -> Generator[Any, None, None]:  # type: ignore
        """Iterate over the values of the Frame."""
        for key in self.keys():
            yield self[key]

    @overload
    def get(self, key: FrameKey[TFrom, TTo]) -> TTo:
        ...

    @overload
    def get(
        self, key: FrameKey[TFrom, TTo], value: Optional[TDefault] = None
    ) -> Union[TTo, TDefault]:
        ...

    @overload
    def get(self, key: str, value: Optional[TDefault] = None) -> Any:
        ...

    def get(self, key: Union[str, FrameKey], value: Optional[TDefault] = None) -> Any:
        """
        Get the value for a given key, returning value if the key is not present.

        :param key: Key to lookup.
        :param value: Default value to use if the key is absent.
        :return: The value for the given key if it is present, otherwise the default
                 return value.
        """
        try:
            return self[key]
        except KeyError:
            return value

    def __ne__(self, other: Any) -> bool:
        return self.raw != other.raw  # type:ignore

    def __repr__(self) -> str:
        rep = "<NarupaFrame"

        for key, value in self.items():
            rep += f" {key}={_print_value(value)}"

        rep += ">"

        return rep


def _print_value(value: Any) -> str:
    if isinstance(value, np.ndarray) and len(value) > 3:
        with np.printoptions(precision=3, suppress=True):
            return f"[{value[0]}, {value[1]}, ... ({len(value)} item(s)]"
    return repr(value)


Mapping.register(NarupaFrame)
