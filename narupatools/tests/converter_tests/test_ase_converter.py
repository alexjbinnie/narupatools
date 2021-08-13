from ase import Atoms

from narupatools.frame import (
    BoxVectors,
    KineticEnergy,
    ParticleCharges,
    ParticleCount,
    ParticleElements,
    ParticleForces,
    ParticleMasses,
    ParticleNames,
    ParticlePositions,
    ParticleResidues,
    ParticleVelocities,
    PotentialEnergy,
    ResidueCount,
    ResidueNames,
)

from .base import BaseTestConverter


class TestASEAtomsConverter(BaseTestConverter):
    frame_to_object_fields = {
        BoxVectors,
        ParticleCharges,
        ParticleCount,
        ParticleElements,
        ParticleForces,
        ParticleMasses,
        ParticlePositions,
        ParticleVelocities,
        PotentialEnergy,
    }

    object_to_frame_fields = {
        BoxVectors,
        KineticEnergy,
        ResidueCount,
        ResidueNames,
        ParticleCharges,
        ParticleCount,
        ParticleElements,
        ParticleForces,
        ParticleMasses,
        ParticleNames,
        ParticlePositions,
        ParticleResidues,
        ParticleVelocities,
        PotentialEnergy,
    }

    @property
    def object_type(self):
        return Atoms
