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
