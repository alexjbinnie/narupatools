from typing import Optional

from .ElementBase import ElementBase

def Element(_tag: str) -> ElementBase: ...
def SubElement(_parent: ElementBase, _tag: str) -> ElementBase: ...
def fromstring(text: str) -> ElementBase: ...
def tostring(
    element_or_tree: ElementBase,
    encoding: Optional[str] = ...,
    method: str = ...,
    pretty_print: bool = ...,
) -> str: ...

__all__ = ["ElementBase", "Element", "fromstring", "tostring"]
