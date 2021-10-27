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

from typing import (
    Any,
    Generic,
    Iterable,
    Iterator,
    Optional,
    Sequence,
    Sized,
    SupportsAbs,
    TypeVar,
)

from .unit import Unit

_T_co = TypeVar("_T_co", covariant=True)
_U = TypeVar("_U")
_TAbs = TypeVar("_TAbs", bound=SupportsAbs)

class Quantity(Generic[_T_co]):
    _value: _T_co
    def __init__(self, value: Optional[_T_co] = ..., unit: Optional[Unit] = ...): ...
    def value_in_unit(self, unit: Unit) -> _T_co: ...
    def __len__(self: Quantity[Sized]) -> int: ...
    def __getitem__(self: Quantity[_U], key: Any) -> Quantity[_U]: ...
    def __iter__(self: Quantity[Iterable[_U]]) -> Iterator[Quantity[_U]]: ...
    def __abs__(self: Quantity[_TAbs]) -> Quantity[_TAbs]: ...
