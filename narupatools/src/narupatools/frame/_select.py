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
        fields |= {ParticleResidues, ResidueNames, ParticleCount, ResidueCount}
    if {"backbone", "nucleicbackbone", "nucleicsugar", "name", "atom"} & keywords:
        fields |= {ParticleNames, ParticleCount}
    if {"type"} & keywords:
        fields |= {ParticleTypes, ParticleCount}
    if {
        "around",
        "sphzone",
        "sphlayer",
        "cyzone",
        "cylayer",
        "point",
        "prop",
    } & keywords:
        fields |= {ParticlePositions, ParticleCount}
    if {"bonded"} & keywords:
        fields |= {BondPairs, BondCount}
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
