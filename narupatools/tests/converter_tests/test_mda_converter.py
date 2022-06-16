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

from MDAnalysis import Universe

from narupatools.frame.fields import (
    BondCount,
    BondPairs,
    BoxVectors,
    ChainCount,
    ChainNames,
    ParticleCharges,
    ParticleCount,
    ParticleElements,
    ParticleForces,
    ParticleMasses,
    ParticleNames,
    ParticlePositions,
    ParticleResidues,
    ParticleTypes,
    ParticleVelocities,
    ResidueChains,
    ResidueCount,
    ResidueIds,
    ResidueNames,
)

from .base import BaseTestConverter


class TestMDAUniverseConverter(BaseTestConverter):
    @property
    def object_type(self):
        return Universe

    frame_to_object_fields = {
        BondCount,
        BondPairs,
        BoxVectors,
        ParticleCharges,
        ParticleCount,
        ParticleElements,
        ParticleForces,
        ParticleMasses,
        ParticleNames,
        ParticlePositions,
        ParticleResidues,
        ParticleVelocities,
        ResidueChains,
        ResidueCount,
        ResidueNames,
    }

    object_to_frame_fields = {
        BondCount,
        BondPairs,
        BoxVectors,
        ChainCount,
        ChainNames,
        ResidueChains,
        ResidueCount,
        ResidueIds,
        ResidueNames,
        ParticleCharges,
        ParticleCount,
        ParticleElements,
        ParticleForces,
        ParticleMasses,
        ParticleNames,
        ParticlePositions,
        ParticleResidues,
        ParticleTypes,
        ParticleVelocities,
    }
