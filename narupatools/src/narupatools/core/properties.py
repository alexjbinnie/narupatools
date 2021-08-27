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

"""Generate properties without having to write getter and setter."""

from typing import Any, Callable, TypeVar

import numpy as np
import numpy.typing as npt

_TValue = TypeVar("_TValue")


def to_property(converter: Callable[[Any], _TValue]) -> Callable[[Callable], _TValue]:
    def _to_property(f: Callable) -> _TValue:
        public_name = f.__name__
        private_name = "_value_" + f.__name__

        def _get(self: Any) -> _TValue:
            try:
                return getattr(self, private_name)  # type: ignore[no-any-return]
            except AttributeError:
                raise AttributeError(f"Property {public_name} of {self} has no value.")

        def _set(self: Any, value: Any) -> None:
            value = converter(value)
            setattr(self, private_name, value)

        def _del(self: Any) -> None:
            delattr(self, private_name)

        return property(fget=_get, fset=_set, fdel=_del)  # type: ignore[return-value]

    return _to_property


def numpy_property(dtype: npt.DTypeLike) -> Callable[[Any], np.ndarray]:
    """Property internally stored as a NumPy array."""

    def to_array(value: Any) -> np.ndarray:
        return np.asarray(value, dtype=dtype)

    return to_property(to_array)


def to_str(value: Any) -> str:
    """Property stored as a string."""
    return str(value)


str_property = to_property(to_str)


def to_float(value: Any) -> float:
    """Property stored as a float."""
    return float(value)


float_property = to_property(to_float)


__all__ = ["numpy_property", "str_property", "float_property"]
