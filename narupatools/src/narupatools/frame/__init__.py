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

"""Code for handling FrameData and other related objects."""

from .converter import convert
from .fields import (
    BondCount,
    BondOrders,
    BondPairs,
    BoxVectors,
    ChainCount,
    ChainNames,
    KineticEnergy,
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
    PotentialEnergy,
    ResidueChains,
    ResidueCount,
    ResidueIds,
    ResidueNames,
    SimulationElapsedSteps,
    SimulationElapsedTime,
    SimulationTotalSteps,
    SimulationTotalTime,
)
from .frame import NarupaFrame
from .patch import *  # noqa: F401, F403
from .trajectory_playback import TrajectoryPlayback

__all__ = [
    "convert",
    "ParticlePositions",
    "ParticleCount",
    "ParticleElements",
    "ParticleNames",
    "ParticleTypes",
    "ParticleResidues",
    "ParticleMasses",
    "ParticleVelocities",
    "ParticleForces",
    "ParticleCharges",
    "ResidueNames",
    "ResidueChains",
    "ResidueCount",
    "ResidueIds",
    "ChainCount",
    "ChainNames",
    "BondCount",
    "BondPairs",
    "BoxVectors",
    "BondOrders",
    "PotentialEnergy",
    "KineticEnergy",
    "SimulationTotalTime",
    "SimulationElapsedTime",
    "SimulationTotalSteps",
    "SimulationElapsedSteps",
    "NarupaFrame",
    "TrajectoryPlayback",
]
