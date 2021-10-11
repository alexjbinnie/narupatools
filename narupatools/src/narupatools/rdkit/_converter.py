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

"""Conversion methods between RDKit and Narupa."""

from typing import Any, Dict, List, Optional, Type, TypeVar, Union

from infinite_sets import InfiniteSet
from MDAnalysis.topology.tables import SYMB2Z, Z2SYMB
from narupa.trajectory import FrameData
from rdkit import Chem

from narupatools.frame import FrameConverter
from narupatools.frame.fields import (
    BondCount,
    BondPairs,
    BondTypes,
    ParticleCount,
    ParticleElements,
    ParticleMasses,
    ParticleNames,
    ParticlePositions,
    ParticleResidues,
    ResidueNames,
)
from narupatools.override import override
from narupatools.physics.units import UnitsNarupa
from narupatools.rdkit._units import UnitsRDKit
from narupatools.util.collections import infinite_seq

_RDKitToNarupa = UnitsRDKit >> UnitsNarupa
_NarupaToRDKit = UnitsNarupa >> UnitsRDKit

ALL_RDKIT_PROPERTIES = frozenset(
    (
        ParticleElements.key,
        ParticlePositions.key,
        ParticleNames.key,
        ParticleMasses.key,
        ParticleCount.key,
        BondPairs.key,
        BondCount.key,
        BondTypes.key,
    )
)

RDKIT_PROPERTIES = frozenset(
    {
        ParticlePositions.key,
        ParticleMasses.key,
        ParticleNames.key,
        ParticleResidues.key,
        ResidueNames.key,
        ParticleElements.key,
        BondPairs.key,
        BondTypes.key,
    }
)

_TType = TypeVar("_TType")

_RDKIT_BOND_TYPES = {
    "a": Chem.BondType.AROMATIC,
    "s": Chem.BondType.SINGLE,
    "d": Chem.BondType.DOUBLE,
    "t": Chem.BondType.TRIPLE,
}

_NT_BOND_TYPES = {
    Chem.BondType.AROMATIC: "a",
    Chem.BondType.SINGLE: "s",
    Chem.BondType.DOUBLE: "d",
    Chem.BondType.TRIPLE: "t",
}


class RDKitConverter(FrameConverter):
    """Frame converter for the RDKit package."""

    @classmethod
    @override(FrameConverter.convert_from_frame)
    def convert_from_frame(  # noqa: D102
        cls,
        frame: FrameData,
        destination: Union[Type[_TType], _TType],
        *,
        fields: InfiniteSet[str],
    ) -> _TType:
        if destination == Chem.rdchem.Mol:
            return frame_to_rdkit_mol(frame)  # type: ignore
        raise NotImplementedError

    @classmethod
    @override(FrameConverter.convert_to_frame)
    def convert_to_frame(  # noqa: D102
        cls,
        object_: _TType,
        /,
        *,
        fields: InfiniteSet[str],
        existing: Optional[FrameData],
    ) -> FrameData:
        if isinstance(object_, Chem.rdchem.Mol):
            return rdkit_mol_to_frame(object_, fields=fields, frame=existing)
        raise NotImplementedError


def rdkit_mol_to_frame(
    mol: Chem.rdchem.Mol,
    fields: InfiniteSet[str] = RDKIT_PROPERTIES,
    frame: Optional[FrameData] = None,
) -> FrameData:
    """
    Convert an RDKit mol to a Narupa FrameData.

    :param mol: RDKit mol to convert.
    :param fields: Fields to populate of the frame data.
    :param frame: Optional preexisting frame data.
    :return: FrameData containing data from RDKit mol.
    """
    if frame is None:
        frame = FrameData()

    indices = {}

    for index, atom in enumerate(mol.GetAtoms()):
        indices[atom.GetIdx()] = index

    if ParticlePositions.key in fields:
        ParticlePositions.set(
            frame, mol.GetConformer().GetPositions() * _RDKitToNarupa.length
        )

    if ParticleMasses.key in fields:
        frame[ParticleMasses] = [atom.GetMass() for atom in mol.GetAtoms()]

    if ParticleElements.key in fields:
        ParticleElements.set(
            frame, [SYMB2Z[atom.GetSymbol()] for atom in mol.GetAtoms()]
        )

    atom_extra = [atom.GetMonomerInfo() for atom in mol.GetAtoms()]

    rdkit_res_number_to_index: Dict[int, int] = {}

    if atom_extra[0] is not None and ParticleNames.key in fields:
        frame[ParticleNames] = [atom.GetName() for atom in atom_extra]

    if atom_extra[0] is not None and (
        ParticleResidues.key in fields or ResidueNames.key in fields
    ):
        res_indices: List[int] = []
        res_names: List[str] = []
        for atom in atom_extra:
            resnum = atom.GetResidueNumber()
            if resnum in rdkit_res_number_to_index:
                res_indices.append(rdkit_res_number_to_index[resnum])
            else:
                res_index = len(res_names)
                rdkit_res_number_to_index[resnum] = res_index
                res_indices.append(res_index)
                res_names.append(atom.GetResidueName())
        if ParticleResidues.key in fields:
            frame[ParticleResidues] = res_indices
        if ResidueNames.key in fields:
            frame[ResidueNames] = res_names

    if BondPairs.key in fields or BondTypes.key in fields:
        bonds = []
        types = []
        for bond in mol.GetBonds():
            id1 = bond.GetBeginAtomIdx()
            id2 = bond.GetEndAtomIdx()
            bonds.append([indices[id1], indices[id2]])
            types.append(_NT_BOND_TYPES.get(bond.GetBondType(), Chem.BondType.SINGLE))
        if BondPairs.key in fields:
            frame[BondPairs] = bonds
        if BondTypes.key in fields:
            frame[BondTypes] = types

    return frame


def frame_to_rdkit_mol(frame: FrameData) -> Chem.rdchem.Mol:
    """
    Convert a Narupa FrameData to an RDKit molecule.

    :param frame: FrameData to convert.
    :return: RDKit molecule.
    """
    mol = Chem.RWMol()

    coordinates = ParticlePositions.get(frame) * _NarupaToRDKit.length
    elements = ParticleElements.get(frame)

    count = len(coordinates)

    ids = []

    for i in range(count):
        atom = Chem.Atom(Z2SYMB[elements[i]])
        idx = mol.AddAtom(atom)
        ids.append(idx)

    if BondPairs.key in frame:
        bonds = BondPairs.get(frame)
        if BondTypes.key in frame:
            types: Any = BondTypes.get(frame)
        else:
            types = infinite_seq("s")
        for i, bond in enumerate(bonds):
            mol.AddBond(
                ids[bond[0]],
                ids[bond[1]],
                _RDKIT_BOND_TYPES.get(types[i], "s"),
            )

    conf = Chem.Conformer(count)
    for i in range(count):
        conf.SetAtomPosition(ids[i], coordinates[i])
    mol.AddConformer(conf)

    mol.UpdatePropertyCache(strict=False)

    Chem.SanitizeMol(mol)

    return mol
