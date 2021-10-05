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
from narupatools.state.typing import Serializable

from ._subgraph import SubgraphObject
from ._typing import Element, Gradient, SingleColor


class Sequence(SubgraphObject):
    @classmethod
    def protein_alpha_carbons(
        cls, **kwargs: Serializable
    ) -> SequenceByPolypeptideAlphaCarbons:
        """
        Generate sequences through contiguous alpha carbons of polypeptides.

        """
        return SequenceByPolypeptideAlphaCarbons(**kwargs)


class SequenceByPolypeptideAlphaCarbons(Sequence):
    """Generate sequences through contiguous alpha carbons of polypeptides."""

    _subgraph_ids: ClassVar[List[str]] = ["polypeptide", "protein"]

    def __init__(self, **kwargs: Serializable):
        """
        Generate sequences through contiguous alpha carbons of polypeptides.

        """
        super().__init__(**kwargs)
