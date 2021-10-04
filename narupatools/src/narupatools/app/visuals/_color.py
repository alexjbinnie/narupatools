from __future__ import annotations

import re
from abc import ABC
from typing import List, Union, Mapping

from narupatools.state import SerializableObject
from narupatools.state.typing import Serializable

_HEX_PATTERN = re.compile("^(?:#|0x)?([A-Fa-f0-9][A-Fa-f0-9])([A-Fa-f0-9][A-Fa-f0-9])([A-Fa-f0-9][A-Fa-f0-9])$")


class Color(SerializableObject, ABC):

    @classmethod
    def from_hex(cls, hex: str):
        return HexColor(hex)

    @classmethod
    def from_named(cls, name: str):
        return NamedColor(name)

    @classmethod
    def from_rgb(cls, r: int, g: int, b: int):
        return ComponentColor([r, g, b])

    @classmethod
    def from_rgba(cls, r: int, g: int, b: int, a: int):
        return ComponentColor([r, g, b, a])

    @classmethod
    def deserialize(cls, value: Serializable) -> Color:
        if isinstance(value, str):
            if _HEX_PATTERN.match(value):
                return cls.from_hex(value)
            return cls.from_named(value)
        elif isinstance(value, List):
            if len(value) in (3, 4):
                return cls.from_rgb(*value)

        raise ValueError(value)

Element = Union[str, int]


class ColorByElement:

    scheme: Union[str, Mapping[Element, Color]]
    """
    Color scheme to be used.
    """

    particle_elements: List[Element]
    """
    List of elements to be used.
    """



class SingleColor(Color):
    pass


class HexColor(SingleColor):

    def __init__(self, hex: str):
        self.hex = hex

    def serialize(self) -> Serializable:
        return self.hex


class NamedColor(SingleColor):

    def __init__(self, name: str):
        self.name = name

    def serialize(self) -> Serializable:
        return self.name


class ComponentColor(SingleColor):

    def __init__(self, components: List[int]):
        self.components = components

    def serialize(self) -> Serializable:
        return self.components
