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
    def copy_file(self, dstfilename: str, overwrite:bool=...): ...
    def close(self) -> None: ...
    def flush(self) -> None: ...
