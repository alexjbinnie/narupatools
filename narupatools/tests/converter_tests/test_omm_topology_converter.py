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

from typing import Type

from openmm.app import Topology

from narupatools.frame.fields import (
    BondCount,
    BondPairs,
    BoxVectors,
    ChainCount,
    ChainNames,
    ParticleCount,
    ParticleElements,
    ParticleNames,
    ParticleResidues,
    ResidueChains,
    ResidueCount,
    ResidueIds,
    ResidueNames,
)

from .base import BaseTestConverter


class TestOpenMMTopologyConverter(BaseTestConverter):
    object_to_frame_fields = {
        BondCount,
        BondPairs,
        BoxVectors,
        ChainCount,
        ChainNames,
        ParticleCount,
        ParticleElements,
        ParticleNames,
        ParticleResidues,
        ResidueChains,
        ResidueCount,
        ResidueIds,
        ResidueNames,
    }

    frame_to_object_fields = {
        BondCount,
        BondPairs,
        BoxVectors,
        ChainCount,
        ChainNames,
        ParticleCount,
        ParticleElements,
        ParticleNames,
        ParticleResidues,
        ResidueChains,
        ResidueCount,
        ResidueIds,
        ResidueNames,
    }

    minimal_fields = {ParticleElements}

    @property
    def object_type(self) -> Type:
        return Topology
