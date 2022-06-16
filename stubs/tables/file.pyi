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
from typing import Any, Literal, Optional, Tuple, Union

from tables.array import Array
from tables.atom import Atom
from tables.earray import EArray
from tables.filters import Filters
from tables.group import Group, RootGroup

def open_file(
    filename: str,
    mode: Literal["r", "w", "a", "r+"] = ...,
    title: str = ...,
    **kwargs: Any
) -> File: ...

class File:
    root: RootGroup
    def create_earray(
        self,
        where: Union[str, Group],
        name: str,
        atom: Optional[Atom] = ...,
        shape: Optional[Tuple[int, ...]] = ...,
        title: str = ...,
        filters: Optional[Filters] = ...,
    ) -> EArray: ...
    def create_group(
        self,
        where: Union[str, Group],
        name: str,
        title: str = ...,
        filters: Optional[Filters] = ...,
    ) -> Group: ...
    def create_array(
        self, where: Union[str, Group], name: str, obj: Any = ..., title: str = ...
    ) -> Array: ...
    def copy_file(self, dstfilename: str, overwrite: bool = ...): ...
    def close(self) -> None: ...
    def flush(self) -> None: ...
