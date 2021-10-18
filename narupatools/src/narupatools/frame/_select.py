from typing import Any

from narupa.trajectory import FrameData

from ._converter import convert
from narupatools.frame.fields import (
    ParticleResidues,
    ResidueNames,
    ParticleNames,
    ResidueIds,
    ParticleTypes,
    ParticlePositions,
    BondPairs, ParticleCount, ResidueCount, BondCount,
)
from MDAnalysis import Universe
import re


def select(frame: Any, /, selection: str):
    keywords = set(re.findall(r"\w+", selection))
    fields = set()
    if {
        "protein",
        "backbone",
        "nucleic",
        "nucleicbackbone",
        "nucleicsugar",
        "resname",
        "atom",
    } & keywords:
        fields |= {ParticleResidues, ResidueNames}
    if {"backbone", "nucleicbackbone", "nucleicsugar", "name", "atom"} & keywords:
        fields |= {ParticleNames}
    if {"type"} & keywords:
        fields |= {ParticleTypes}
    if {
        "around",
        "sphzone",
        "sphlayer",
        "cyzone",
        "cylayer",
        "point",
        "prop",
    } & keywords:
        fields |= {ParticlePositions}
    if {"bonded"} & keywords:
        fields |= {BondPairs}

    if not isinstance(frame, FrameData):
        frame = convert(frame, FrameData, fields=fields)

    if missing_keys := fields - set(frame.keys()):
        raise ValueError(
            f"Cannot use selection '{selection}' - key(s) {missing_keys} are unavailable."
        )
    universe = convert(frame, Universe, fields=fields)
    return universe.select_atoms(selection).ix


FrameData.select = select
