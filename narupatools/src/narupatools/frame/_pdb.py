from io import StringIO
from typing import Any, Dict, List, Optional

import numpy as np
from MDAnalysis.topology.tables import Z2SYMB
from narupa.trajectory import FrameData

from narupatools.frame.fields import (
    BondPairs,
    BoxVectors,
    ParticleCharges,
    ParticleElements,
    ParticleNames,
    ParticlePositions,
    ParticleResidues,
    ResidueNames,
)
from narupatools.physics.units import angstrom, degree, meter, nano, radian
from narupatools.physics.vector import angle, magnitude
from narupatools.util.collections import infinite_seq

PDB_STD_RESIDUES = frozenset(
    {
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "PYL",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "SEC",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
        "ASX",
        "GLX",
        "UNK",
        "A",
        "C",
        "G",
        "I",
        "U",
        "DA",
        "DC",
        "DG",
        "DI",
        "DT",
        "DU",
        "N",
    }
)
"""
Set of standard residues as defined by the PDB file format.

Any atom in a residue with one of these names should be stored under an ATOM record
instead of a HETATM record.
"""


def frame_to_pdb_string(frame_data: FrameData) -> str:
    """
    Convert a Narupa FrameData to the contents of a PDB file.

    This is useful for interacting with processes which accept a PDB file as input.

    :param frame_data: FrameData to convert.
    :return: Contents of a PDB file as a string.
    """
    file = StringIO()

    try:
        positions: Any = ParticlePositions.get(frame_data) * 10.0
    except KeyError:
        positions = []
    count = len(positions)
    atomnames = ParticleNames.get_with_default(frame_data, infinite_seq(""))
    elements = ParticleElements.get_with_default(frame_data, infinite_seq(0))
    atomresidues = ParticleResidues.get_with_default(frame_data, infinite_seq(0))
    resnames = ResidueNames.get_with_default(frame_data, infinite_seq("UNL"))
    try:
        charges = ParticleCharges.get(frame_data)
    except KeyError:
        charges = np.zeros(count)

    remark = "REMARK GENERATED BY NARUPATOOLS"
    file.write(f"{remark:<80}\n")

    # Write CRYST1 entry to store simulation box.
    try:
        box_vectors: Optional[np.ndarray] = BoxVectors.get(frame_data)
    except KeyError:
        box_vectors = None

    if box_vectors is not None:
        a = magnitude(box_vectors[0]) * (nano * meter >> angstrom)
        b = magnitude(box_vectors[1]) * (nano * meter >> angstrom)
        c = magnitude(box_vectors[2]) * (nano * meter >> angstrom)
        if not np.isclose(a, 0) and not np.isclose(b, 0) and not np.isclose(c, 0):
            alpha = angle(box_vectors[1], box_vectors[2]) * (radian >> degree)
            beta = angle(box_vectors[0], box_vectors[2]) * (radian >> degree)
            gamma = angle(box_vectors[0], box_vectors[1]) * (radian >> degree)
            sgroup = "P 1"
            file.write(
                f"CRYST1{a:9.3f}{b:9.3f}{c:9.3f}{alpha:7.2f}{beta:7.2f}{gamma:7.2f} {sgroup:<11}   1          \n"
            )

    for index, name, residue, position, element, charge in zip(
        range(count), atomnames, atomresidues, positions, elements, charges
    ):
        record_type = "ATOM  " if resnames[residue] in PDB_STD_RESIDUES else "HETATM"
        symb = Z2SYMB[element] if element > 0 else "X"
        file.write(
            f"{record_type}{index + 1:5} {name:4.4} {resnames[residue]:3.3} "
            f"A{residue + 1:4}    {position[0]:8.3f}{position[1]:8.3f}"
            f"{position[2]:8.3f}  1.00  0.00          {symb:>2}"
            f"{int(charge):2}\n"
        )
    bonded: Dict[int, List[int]] = {}
    for bond in BondPairs.get_with_default(frame_data, []):
        if bond[0] not in bonded:
            bonded[bond[0]] = []
        if bond[1] not in bonded:
            bonded[bond[1]] = []
        bonded[bond[0]].append(bond[1])
        bonded[bond[1]].append(bond[0])

    for index in range(count):
        if index not in bonded:
            continue
        line = f"CONECT{index + 1:5}"
        for n in bonded[index]:
            line += f"{n + 1:5}"
        file.write(f"{line:<80}\n")

    return file.getvalue()