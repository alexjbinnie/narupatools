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

"""Conversion methods between MDAnalysis and Narupa."""

import contextlib
from typing import List, Optional, Type, TypeVar, Union

import numpy as np
from MDAnalysis import AtomGroup, NoDataError, Universe
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.core.topology import Topology
from MDAnalysis.core.topologyattrs import (
    Atomnames,
    Bonds,
    Charges,
    Elements,
    Masses,
    Resnames,
    TopologyAttr,
)
from ase.geometry import cell_to_cellpar, cellpar_to_cell
from infinite_sets import InfiniteSet, everything
from narupa.trajectory import FrameData

from narupatools.core.units import UnitsNarupa
from narupatools.frame import (
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
from narupatools.frame.converter import FrameConverter
from narupatools.frame.utils import atomic_numbers_to_symbols
from narupatools.mdanalysis.units import UnitsMDAnalysis
from narupatools.mdanalysis.utils import guess_atomic_number

MDAnalysisToNarupa = UnitsMDAnalysis >> UnitsNarupa
NarupaToMDAnalysis = UnitsNarupa >> UnitsMDAnalysis

ALL_MDA_PROPERTIES = frozenset(
    (
        ParticleCharges.key,
        ParticleElements.key,
        ParticlePositions.key,
        ParticleForces.key,
        ParticleVelocities.key,
        ParticleNames.key,
        ParticleMasses.key,
        ParticleCount.key,
        ResidueCount.key,
        ResidueNames.key,
        ResidueChains.key,
        ResidueIds.key,
        ParticleTypes.key,
        ParticleResidues.key,
        ParticleTypes.key,
        ChainCount.key,
        ChainNames.key,
        BondPairs.key,
        BoxVectors.key,
        BondCount.key,
    )
)

MDA_PROPERTIES = frozenset(
    (
        ParticleElements.key,
        ParticlePositions.key,
        ParticleNames.key,
        ParticleCount.key,
        ResidueCount.key,
        ResidueNames.key,
        ResidueChains.key,
        ResidueIds.key,
        ParticleTypes.key,
        ParticleResidues.key,
        ParticleTypes.key,
        ChainCount.key,
        ChainNames.key,
        BondPairs.key,
        BoxVectors.key,
        BondCount.key,
    )
)

_TType = TypeVar("_TType")


class MDAnalysisConverter(FrameConverter):
    """Frame converter for the MDAnalysis package."""

    @classmethod
    def convert_from_frame(  # noqa: D102
        cls,
        frame: FrameData,
        type: Union[Type[_TType], _TType],
        *,
        fields: InfiniteSet[str],
    ) -> _TType:
        if type == Topology:
            return frame_to_mdanalysis_topology(frame, fields=fields)  # type: ignore
        if type == Universe:
            return frame_to_mdanalysis_universe(frame, fields=fields)  # type: ignore
        raise NotImplementedError()

    @classmethod
    def convert_to_frame(  # noqa: D102
        cls, object: _TType, *, fields: InfiniteSet[str], existing: Optional[FrameData]
    ) -> FrameData:
        if isinstance(object, Universe):
            return mdanalysis_universe_to_frame(object, fields=fields, frame=existing)
        if isinstance(object, AtomGroup):
            return mdanalysis_atomgroup_to_frame(object, fields=fields, frame=existing)
        raise NotImplementedError()


ALL_FIELDS = everything()


def frame_to_mdanalysis_universe(
    frame: FrameData, *, fields: InfiniteSet[str] = ALL_FIELDS
) -> Universe:
    """
    Convert a Narupa frame to an MDAnlysis universe.

    :param frame: Narupa frame to convert.
    :param fields: Fields to copy from the frame data.
    :return: Universe with the given data included.
    """
    topology = frame_to_mdanalysis_topology(frame, fields=fields)
    universe = Universe(topology)

    coords = np.zeros((0, 3))
    velocities = None
    forces = None
    dimensions = None

    if ParticlePositions.key in fields:
        with contextlib.suppress(KeyError):
            coords = ParticlePositions.get(frame) * NarupaToMDAnalysis.length

    if ParticleVelocities.key in fields:
        with contextlib.suppress(KeyError):
            velocities = ParticleVelocities.get(frame) * NarupaToMDAnalysis.velocity

    if ParticleForces.key in fields:
        with contextlib.suppress(KeyError):
            forces = ParticleForces.get(frame) * NarupaToMDAnalysis.force

    if BoxVectors.key in fields:
        with contextlib.suppress(KeyError):
            dimensions = cell_to_cellpar(
                BoxVectors.get(frame) * NarupaToMDAnalysis.length
            )

    universe.trajectory = MemoryReader(
        coords, order="fac", velocities=velocities, forces=forces, dimensions=dimensions
    )

    return universe


def frame_to_mdanalysis_topology(  # noqa: C901
    frame: FrameData, *, fields: InfiniteSet[str] = ALL_FIELDS
) -> Topology:
    """
    Convert a Narupa frame to an MDAnalysis topology.

    :param frame: Narupa frame to convert.
    :param fields: Set of fields to copy.
    :return: MDAnalysis topology with the given fields copied over.
    """
    try:
        natoms = ParticleCount.get(frame)
    except KeyError:
        natoms = 0

    try:
        nress = ResidueCount.get(frame)
    except KeyError:
        nress = 0

    try:
        nsegs = ChainCount.get(frame)
    except KeyError:
        nsegs = 1

    attrs: List[TopologyAttr] = []

    if ResidueNames.key in fields:
        with contextlib.suppress(KeyError):
            attrs.append(Resnames(ResidueNames.get(frame)))

    if ParticleNames.key in fields:
        with contextlib.suppress(KeyError):
            attrs.append(Atomnames(ParticleNames.get(frame)))

    if ParticleCharges.key in fields:
        with contextlib.suppress(KeyError):
            attrs.append(Charges(ParticleCharges.get(frame)))

    if ParticleMasses.key in fields:
        with contextlib.suppress(KeyError):
            attrs.append(Masses(ParticleMasses.get(frame)))

    if ParticleElements.key in fields:
        with contextlib.suppress(KeyError):
            atomatomicnumbers = ParticleElements.get(frame)
            atomelements = atomic_numbers_to_symbols(atomatomicnumbers)
            attrs.append(Elements(atomelements))

    res_segidx = np.array([0] * nress)
    if ResidueChains.key in fields:
        with contextlib.suppress(KeyError):
            res_segidx = ResidueChains.get(frame)

    atom_residx = np.array([0] * natoms)
    if ParticleResidues.key in fields:
        with contextlib.suppress(KeyError):
            atom_residx = ParticleResidues.get(frame)

    if BondPairs.key in fields:
        with contextlib.suppress(KeyError):
            attrs.append(Bonds(BondPairs.get(frame)))

    return Topology(
        n_atoms=natoms,
        n_res=nress,
        n_seg=nsegs,
        attrs=attrs,
        residue_segindex=res_segidx,
        atom_resindex=atom_residx,
    )


def mdanalysis_universe_to_frame(
    universe: Universe,
    *,
    fields: InfiniteSet[str] = MDA_PROPERTIES,
    frame: Optional[FrameData] = None,
) -> FrameData:
    """
    Convert an MDAnalysis Universe into a Narupa FrameData.

    :param universe: The MDAnalysis universe to convert.
    :param fields: Fields to add to FrameData.
    :param frame: Optional pre-existing FrameData to populate.
    :return: FrameData populated with requested fields.
    """
    return mdanalysis_atomgroup_to_frame(universe.atoms, fields=fields, frame=frame)


def mdanalysis_atomgroup_to_frame(
    group: AtomGroup,
    *,
    fields: InfiniteSet[str] = MDA_PROPERTIES,
    frame: Optional[FrameData] = None,
) -> FrameData:
    """
    Convert an MDAnalysis AtomGroup into a Narupa FrameData.

    :param group: The MDAnalysis AtomGroup to convert.
    :param fields: Fields to add to FrameData.
    :param frame: Optional pre-existing FrameData to populate.
    :return: FrameData populated with requested fields.
    """
    if frame is None:
        frame = FrameData()
    _add_mda_atoms_core_attrs(group, fields, frame)
    _add_mda_atoms_extra_attrs(group, fields, frame)
    _add_mda_residue_attrs(group, fields, frame)
    _add_mda_segments_attrs(group, fields, frame)
    _add_mda_bonds_attrs(group, fields, frame)
    _add_mda_box(group, fields, frame)
    return frame


def _is_atom_group_selection(group: AtomGroup) -> bool:
    return len(group) != len(group.universe.atoms)


def _add_mda_atoms_core_attrs(
    group: AtomGroup, fields: InfiniteSet[str], frame: FrameData
) -> None:
    if ParticleCount.key in fields:
        ParticleCount.set(frame, group.n_atoms)
    if ParticlePositions.key in fields:
        with contextlib.suppress(NoDataError):
            ParticlePositions.set(frame, group.positions * MDAnalysisToNarupa.length)
    if ParticleNames.key in fields:
        with contextlib.suppress(NoDataError):
            ParticleNames.set(frame, group.names)
    if ParticleElements.key in fields:
        with contextlib.suppress(NoDataError):
            ParticleElements.set(frame, guess_atomic_number(group))
    if ParticleTypes.key in fields:
        with contextlib.suppress(NoDataError):
            ParticleTypes.set(frame, group.types)
    if ParticleResidues.key in fields:
        with contextlib.suppress(NoDataError):
            if _is_atom_group_selection(group):
                residue_ix_to_index = {
                    residue_ix: index
                    for index, residue_ix in enumerate(group.residues.ix)
                }
                ParticleResidues.set(
                    frame, [residue_ix_to_index[ix] for ix in group.resindices]
                )
            else:
                ParticleResidues.set(frame, group.resindices)


def _add_mda_atoms_extra_attrs(
    group: AtomGroup, fields: InfiniteSet[str], frame: FrameData
) -> None:
    if ParticleVelocities.key in fields:
        with contextlib.suppress(NoDataError):
            ParticleVelocities.set(
                frame, group.velocities * MDAnalysisToNarupa.velocity
            )
    if ParticleForces.key in fields:
        with contextlib.suppress(NoDataError):
            ParticleForces.set(frame, group.forces * MDAnalysisToNarupa.force)
    if ParticleCharges.key in fields:
        with contextlib.suppress(NoDataError):
            ParticleCharges.set(frame, group.charges * MDAnalysisToNarupa.charge)
    if ParticleMasses.key in fields:
        with contextlib.suppress(NoDataError):
            ParticleMasses.set(frame, group.masses * MDAnalysisToNarupa.mass)


def _add_mda_residue_attrs(
    group: AtomGroup, fields: InfiniteSet[str], frame: FrameData
) -> None:
    if ResidueCount.key in fields:
        ResidueCount.set(frame, group.n_residues)
    if ResidueNames.key in fields:
        with contextlib.suppress(NoDataError):
            ResidueNames.set(frame, group.residues.resnames)
    if ResidueIds.key in fields:
        with contextlib.suppress(NoDataError):
            ResidueIds.set(frame, [str(i) for i in group.residues.resids])
    if ResidueChains.key in fields:
        with contextlib.suppress(NoDataError):
            if _is_atom_group_selection(group):
                segment_ix_to_index = {
                    segment_ix: index
                    for index, segment_ix in enumerate(group.segments.ix)
                }
                ResidueChains.set(
                    frame, [segment_ix_to_index[ix] for ix in group.residues.segindices]
                )
            else:
                ResidueChains.set(frame, group.residues.segindices)


def _add_mda_segments_attrs(
    group: AtomGroup, fields: InfiniteSet[str], frame: FrameData
) -> None:
    if ChainCount.key in fields:
        ChainCount.set(frame, group.n_segments)
    if ChainNames.key in fields:
        with contextlib.suppress(NoDataError):
            ChainNames.set(frame, group.segments.segids)


def _add_mda_bonds_attrs(
    group: AtomGroup, fields: InfiniteSet[str], frame: FrameData
) -> None:
    if _is_atom_group_selection(group) and (
        BondPairs.key in fields or BondCount.key in fields
    ):
        particle_ix_to_index = {
            particle_ix: index for index, particle_ix in enumerate(group.ix)
        }
        bond_pairs = []
        for bond in group.bonds.indices:
            with contextlib.suppress(KeyError):
                bond_pairs.append(
                    [particle_ix_to_index[bond[0]], particle_ix_to_index[bond[1]]]
                )
        if BondPairs.key in fields:
            BondPairs.set(frame, bond_pairs)
        if BondCount.key in fields:
            BondCount.set(frame, len(bond_pairs))
    else:
        if BondCount.key in fields:
            with contextlib.suppress(NoDataError):
                BondCount.set(frame, len(group.bonds))
        if BondPairs.key in fields:
            with contextlib.suppress(NoDataError):
                BondPairs.set(frame, group.bonds.indices)


def _add_mda_box(group: AtomGroup, fields: InfiniteSet[str], frame: FrameData) -> None:
    if BoxVectors.key in fields and np.any(group.dimensions):
        BoxVectors.set(
            frame, cellpar_to_cell(group.dimensions) * MDAnalysisToNarupa.length
        )
