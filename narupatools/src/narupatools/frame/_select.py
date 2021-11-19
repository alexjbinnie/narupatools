import re
from typing import Any, Set

import numpy as np
from MDAnalysis import Universe
from narupa.trajectory import FrameData

from narupatools.frame.fields import (
    BondCount,
    BondPairs,
    ParticleCount,
    ParticleElements,
    ParticleNames,
    ParticlePositions,
    ParticleResidues,
    ParticleTypes,
    ResidueCount,
    ResidueNames,
)

from ._converter import convert
from ._frame_source import FrameSource, TrajectorySource


def _get_selection_fields(selection: str, /) -> Set[str]:
    """Get the set of fields required to determine a selection."""
    keywords = set(re.findall(r"\w+", selection))
    fields: Set[str] = set()
    if {
        "protein",
        "backbone",
        "nucleic",
        "nucleicbackbone",
        "nucleicsugar",
        "resname",
        "atom",
    } & keywords:
        fields |= {ParticleResidues, ResidueNames, ParticleCount, ResidueCount}  # type: ignore
    if {"backbone", "nucleicbackbone", "nucleicsugar", "name", "atom"} & keywords:
        fields |= {ParticleNames, ParticleCount}  # type: ignore
    if {"type"} & keywords:
        fields |= {ParticleTypes, ParticleCount}  # type: ignore
    if {
        "around",
        "sphzone",
        "sphlayer",
        "cyzone",
        "cylayer",
        "point",
        "prop",
    } & keywords:
        fields |= {ParticlePositions, ParticleCount}  # type: ignore
    if {"bonded"} & keywords:
        fields |= {BondPairs, BondCount}  # type: ignore
    if {"element"} & keywords:
        fields |= {ParticleElements}
    return fields


def select(frame: Any, /, selection: str) -> np.ndarray:
    """Calculate the particle indices in a FrameData given by a selection."""
    fields = _get_selection_fields(selection)

    if isinstance(frame, FrameSource):
        frame = frame.get_frame(fields=fields)
    elif isinstance(frame, TrajectorySource):
        frame = frame.get_frame(fields=fields, index=0)
    elif not isinstance(frame, FrameData):
        frame = convert(frame, FrameData, fields=fields)

    if missing_keys := fields - set(frame.keys()):
        raise ValueError(
            f"Cannot use selection '{selection}' - key(s) {missing_keys} are unavailable."
        )
    universe = convert(frame, Universe, fields=fields)
    return universe.select_atoms(selection).ix
