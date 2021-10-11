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

# mypy: ignore-errors

"""Generate properties without having to write getter and setter."""
import collections
import contextlib
import typing
from typing import Any, Callable, Optional, Tuple, TypeVar

import numpy as np
import numpy.typing as npt

from narupatools.physics import quaternion
from narupatools.physics.transformation import Rotation
from narupatools.state import SerializableObject

_TValue = TypeVar("_TValue")


def to_property(converter: Callable[[Any], _TValue]) -> Callable[[Callable], _TValue]:
    def _to_property(f: Callable) -> _TValue:
        public_name = f.__name__
        private_name = "_value_" + f.__name__

        def _get(self: Any) -> _TValue:
            try:
                return getattr(self, private_name)  # type: ignore[no-any-return]
            except AttributeError as e:
                raise AttributeError(
                    f"Property {public_name} of {self} has no value."
                ) from e

        def _set(self: Any, value: Any) -> None:
            value = converter(value)
            setattr(self, private_name, value)

        def _del(self: Any) -> None:
            delattr(self, private_name)

        return property(fget=_get, fset=_set, fdel=_del, doc=f.__doc__)  # type: ignore[return-value]

    return _to_property


def numpy(
    *, dtype: npt.DTypeLike, shape: Tuple[Optional[int]]
) -> Callable[[Any], np.ndarray]:
    """Property internally stored as a NumPy array."""

    def to_array(value: Any) -> np.ndarray:
        v = np.asarray(value, dtype=dtype)
        if len(v.shape) != len(shape):
            raise TypeError(f"Incompatible shapes {v.shape} and {shape}")
        for dim, dim_expected in zip(v.shape, shape):
            if dim_expected and dim != dim_expected:
                raise TypeError(f"Incompatible shapes {v.shape} and {shape}")
        return v

    return to_property(to_array)


def to_str(value: Any) -> str:
    """Property stored as a string."""
    return str(value)


def to_bool(value: Any) -> bool:
    """Property stored as a bool."""
    return bool(value)


string = to_property(to_str)

boolean = to_property(to_bool)


def to_float(value: Any) -> float:
    """Property stored as a float."""
    return float(value)


def to_int(value: Any) -> int:
    """Property stored as an int."""
    return int(value)


def to_unit_quaternion(value: Any) -> quaternion:
    if isinstance(value, quaternion):
        return value
    if isinstance(value, Rotation):
        return value.versor
    if len(value) != 4:
        raise ValueError(f"Can't interperate {value} as quaternion")
    return quaternion(*value)


number = to_property(to_float)

integer = to_property(to_int)

unit_quaternion = to_property(to_unit_quaternion)


def to_list(item_converter):
    def _conv(value):
        return [item_converter(item) for item in value]

    return _conv


def to_serializable(serializable_type):
    def _conv(value):
        if isinstance(value, serializable_type):
            return value
        return serializable_type.deserialize(value)

    return _conv


def get_converter(annot):
    if annot == float:
        return to_float
    if annot == int:
        return to_int
    if annot == str:
        return to_str
    if typing.get_origin(annot) == list:
        return to_list(get_converter(typing.get_args(annot)[0]))
    if typing.get_origin(annot) == collections.abc.Mapping:
        return lambda x: x
    if isinstance(annot, type) and issubclass(annot, SerializableObject):
        return to_serializable(annot)
    if typing.get_origin(annot) == typing.Union:
        arg_types = list(typing.get_args(annot))
        arg_types.sort(
            key=lambda x: 0
            if (isinstance(x, type) and issubclass(x, SerializableObject))
            else 1
        )
        converters = [get_converter(arg) for arg in arg_types]
        arg_types = [arg for arg in arg_types if len(typing.get_args(arg)) == 0]

        def to_union(value):
            for arg in arg_types:
                if isinstance(value, arg):
                    return value
            for conv in converters:
                with contextlib.suppress(ValueError, TypeError):
                    return conv(value)
            raise ValueError

        return to_union
    return lambda x: x


def auto(f: Callable) -> Any:
    """Convert an empty annotated function into a property."""
    public_name = f.__name__
    private_name = "_value_" + f.__name__

    converter = get_converter(typing.get_type_hints(f)["return"])

    def _get(self: Any) -> _TValue:
        try:
            return getattr(self, private_name)  # type: ignore[no-any-return]
        except AttributeError as e:
            raise AttributeError(
                f"Property {public_name} of {self} has no value."
            ) from e

    def _set(self: Any, value: Any) -> None:
        value = converter(value)
        setattr(self, private_name, value)

    def _del(self: Any) -> None:
        delattr(self, private_name)

    return property(fget=_get, fset=_set, fdel=_del, doc=f.__doc__)  # type: ignore[return-value]


__all__ = ["numpy", "string", "boolean", "number", "unit_quaternion", "auto"]
