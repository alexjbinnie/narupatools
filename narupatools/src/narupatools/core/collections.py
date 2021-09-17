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

"""Miscellaneous collections for internal use."""

from __future__ import annotations

from typing import Any, Generic, TypeVar

_T = TypeVar("_T")


class infinite_seq(Generic[_T]):
    """
    Represents an (infinitely long) sequence of the same value.

    This is similar to :class:`itertools.repeat`, except it also returns the value when called using
    index notation.

    .. code-block:: python

        seq = infinite_seq("my_string")

        seq[242]  # returns "my_string"

    """

    def __init__(self, value: _T):
        self._value = value

    def __getitem__(self, i: Any) -> _T:
        return self._value

    def __iter__(self) -> infinite_seq[_T]:
        return self

    def __next__(self) -> _T:
        return self._value
