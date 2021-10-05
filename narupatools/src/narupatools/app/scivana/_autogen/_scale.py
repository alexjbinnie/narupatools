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


class Scale(SubgraphObject):
    @classmethod
    def covalent(
        cls, *, scale: Optional[float] = None, **kwargs: Serializable
    ) -> ScaleByCovalent:
        """
        Scale each atom by its covalent radii.

        :param scale: Scaling factor to apply to each covalent radius.
        """
        return ScaleByCovalent(scale=scale, **kwargs)

    @classmethod
    def by_secondary_stucture(
        cls,
        *,
        scheme: Optional[Union[str, Mapping[str, float]]] = None,
        **kwargs: Serializable
    ) -> ScaleBySecondaryStructure:
        """
        Scale each particle by its secondary structure assignment.

        :param scheme: Scheme mapping secondary structure assignments to scales.
        """
        return ScaleBySecondaryStructure(scheme=scheme, **kwargs)

    @classmethod
    def van_der_waals(
        cls, *, scale: Optional[float] = None, **kwargs: Serializable
    ) -> ScaleByVanDerWaals:
        """
        Scale each atom by its Van der Waals radii.

        :param scale: Scaling factor to apply to each Van der Waals radius.
        """
        return ScaleByVanDerWaals(scale=scale, **kwargs)


class ScaleByCovalent(Scale):
    """Scale each atom by its covalent radii."""

    _subgraph_ids: ClassVar[List[str]] = ["covalent"]

    def __init__(self, *, scale: Optional[float] = None, **kwargs: Serializable):
        """
        Scale each atom by its covalent radii.

        :param scale: Scaling factor to apply to each covalent radius.
        """
        super().__init__(**kwargs)
        if scale is not None:
            self.scale = scale

    @properties.auto
    def scale(self) -> float:
        """Scaling factor to apply to each covalent radius."""


class ScaleBySecondaryStructure(Scale):
    """Scale each particle by its secondary structure assignment."""

    _subgraph_ids: ClassVar[List[str]] = ["secondary structure"]

    def __init__(
        self,
        *,
        scheme: Optional[Union[str, Mapping[str, float]]] = None,
        **kwargs: Serializable
    ):
        """
        Scale each particle by its secondary structure assignment.

        :param scheme: Scheme mapping secondary structure assignments to scales.
        """
        super().__init__(**kwargs)
        if scheme is not None:
            self.scheme = scheme

    @properties.auto
    def scheme(self) -> Union[str, Mapping[str, float]]:
        """Scheme mapping secondary structure assignments to scales."""


class ScaleByVanDerWaals(Scale):
    """Scale each atom by its Van der Waals radii."""

    _subgraph_ids: ClassVar[List[str]] = ["vdw"]

    def __init__(self, *, scale: Optional[float] = None, **kwargs: Serializable):
        """
        Scale each atom by its Van der Waals radii.

        :param scale: Scaling factor to apply to each Van der Waals radius.
        """
        super().__init__(**kwargs)
        if scale is not None:
            self.scale = scale

    @properties.auto
    def scale(self) -> float:
        """Scaling factor to apply to each Van der Waals radius."""
