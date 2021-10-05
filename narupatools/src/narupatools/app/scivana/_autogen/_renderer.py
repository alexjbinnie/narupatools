from __future__ import annotations

from typing import Optional, Union

import narupatools.util.properties as properties
from narupatools.state import SharedStateObject
from narupatools.state.typing import Serializable

from ._color import Color
from ._render import Render
from ._scale import Scale
from ._sequence import Sequence
from ._typing import SingleColor


class Renderer(SharedStateObject):
    @properties.auto
    def color(self) -> Union[Color, SingleColor]:
        """"""

    @properties.auto
    def scale(self) -> Union[float, Scale]:
        """"""

    @properties.auto
    def width(self) -> Union[float, Scale]:
        """"""

    @properties.auto
    def sequence(self) -> Sequence:
        """"""

    @properties.auto
    def render(self) -> Render:
        """"""

    @classmethod
    def create(
        cls,
        color: Optional[Union[Color, SingleColor]] = None,
        scale: Optional[Union[float, Scale]] = None,
        width: Optional[Union[float, Scale]] = None,
        render: Optional[Render] = None,
        sequence: Optional[Sequence] = None,
        **kwargs: Serializable
    ) -> Renderer:
        """
        Create a new custom visualiser.

        :param color: Color to use for the visualiser.
        :param scale: Scale to use for the visualiser.
        :param width: Widths to use for the visualiser if required.
        :param render: Render style to use for the visualiser.
        :param sequence: Sequences to use for the visualiser if required.
        """
        return Renderer(
            color=color,
            scale=scale,
            width=width,
            sequence=sequence,
            render=render,
            **kwargs
        )
