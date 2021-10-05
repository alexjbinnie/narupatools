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

from __future__ import annotations

from typing import ClassVar, List, Mapping, Optional, Union

import narupatools.util.properties as properties
from narupatools.state import SharedStateObject
from narupatools.state.typing import Serializable

from ._color import Color
from ._render import Render
from ._scale import Scale
from ._sequence import Sequence
from ._subgraph import SubgraphObject
from ._typing import Element, Gradient, SingleColor


class Renderer(SharedStateObject):
    @properties.auto
    def color(self) -> Union[Color, SingleColor]:
        """ """

    @properties.auto
    def scale(self) -> Union[float, Scale]:
        """ """

    @properties.auto
    def width(self) -> Union[float, Scale]:
        """ """

    @properties.auto
    def sequence(self) -> Sequence:
        """ """

    @properties.auto
    def render(self) -> Render:
        """ """

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
