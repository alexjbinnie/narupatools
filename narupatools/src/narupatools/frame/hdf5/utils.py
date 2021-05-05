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

"""Utilities for handling the NarupaTools HDF5 Format."""

import itertools
import json
from typing import Any

import numpy as np
from MDAnalysis.topology.tables import Z2SYMB
from narupa.trajectory import FrameData

from narupatools.frame import (
    BondPairs,
    ChainCount,
    ParticleCount,
    ParticleElements,
    ParticleNames,
    ParticleResidues,
    ResidueChains,
    ResidueCount,
    ResidueNames,
)


def generate_topology(frame: FrameData, /) -> bytes:
    """
    Generate ASCII-encoded JSON string of the topology of the provided Frame.

    :param frame: Frame to generate topology for.
    """
    bonds: Any = BondPairs.get_with_default(frame, default=np.empty((0,))).tolist()
    chain_count = ChainCount.get_with_default(frame, calculate=True, default=1)

    residue_count = ResidueCount.get_with_default(frame, default=1, calculate=True)
    residue_names = ResidueNames.get_with_default(frame, default=itertools.repeat(""))
    residue_chains = ResidueChains.get_with_default(frame, default=itertools.repeat(0))

    atom_count = ParticleCount.get(frame, calculate=True)
    atom_names = ParticleNames.get_with_default(frame, default=itertools.repeat(""))
    atom_elements = ParticleElements.get_with_default(
        frame, calculate=True, default=itertools.repeat(None)
    )
    atom_residues = ParticleResidues.get_with_default(
        frame, default=itertools.repeat(0)
    )

    chains = []
    residues = []

    for i in range(chain_count):
        chains.append({"index": i, "residues": []})

    for res_index, res_name, res_chain in zip(
        range(residue_count), residue_names, residue_chains
    ):
        res = {
            "index": res_index,
            "resSeq": res_index + 1,
            "name": res_name,
            "atoms": [],
        }
        residues.append(res)
        chains[res_chain]["residues"].append(res)

    for atom_index, atom_name, atom_element, atom_res in zip(
        range(atom_count), atom_names, atom_elements, atom_residues
    ):
        residues[atom_res]["atoms"].append(
            {
                "index": atom_index,
                "name": atom_name,
                "element": Z2SYMB[atom_element]
                if (atom_element is not None and atom_element != 0)
                else "",
            }
        )

    topology = {"bonds": bonds, "chains": chains}

    return json.dumps(topology).encode("ascii")
