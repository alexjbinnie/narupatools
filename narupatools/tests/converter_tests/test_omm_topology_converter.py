from typing import Type

from simtk.openmm.app import Topology

from narupatools.frame import (
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
