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

"""Conversion methods for converting between MDTraj and Narupa."""

from typing import Optional, Type, TypeVar, Union

from infinite_sets import InfiniteSet
from mdtraj import Topology, Trajectory
from narupa.trajectory.frame_data import FrameData

from narupatools.frame import FrameConverter
from narupatools.frame.fields import (
    BondCount,
    BondPairs,
    BoxVectors,
    ChainCount,
    ParticleCount,
    ParticleElements,
    ParticleNames,
    ParticlePositions,
    ParticleResidues,
    ResidueChains,
    ResidueCount,
    ResidueNames,
)
from narupatools.mdtraj._units import UnitsMDTraj
from narupatools.override import override
from narupatools.physics.units import UnitsNarupa

MDTRAJ_PROPERTIES = frozenset(
    (
        ParticlePositions,
        BondPairs,
        ParticleResidues,
        ParticleElements,
        ParticleNames,
        ResidueNames,
        ResidueChains,
        ParticleCount,
        ResidueCount,
        BondCount,
        ChainCount,
        BoxVectors,
    )
)

_MDTrajToNarupa = UnitsMDTraj >> UnitsNarupa

_TType = TypeVar("_TType")


class MDTrajConverter(FrameConverter):
    """FrameConverter for the mdtraj package."""

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
        if isinstance(object_, Trajectory):
            return mdtraj_trajectory_to_frame(object_, fields=fields, frame=existing)
        if isinstance(object_, Topology):
            return mdtraj_topology_to_frame(object_, fields=fields, frame=existing)
        raise NotImplementedError

    @classmethod
    @override(FrameConverter.convert_from_frame)
    def convert_from_frame(  # noqa: D102
        cls,
        frame: FrameData,
        destination: Union[Type[_TType], _TType],
        *,
        fields: InfiniteSet[str],
    ) -> _TType:
        raise NotImplementedError


def mdtraj_trajectory_to_frame(
    trajectory: Trajectory,
    *,
    frame_index: int = 0,
    fields: InfiniteSet[str] = MDTRAJ_PROPERTIES,
    frame: Optional[FrameData] = None,
) -> FrameData:
    """
    Convert an MDTraj trajectory into a FrameData.

    :param trajectory: Trajectory to convert.
    :param frame_index: Index of the frame in the trajectory to convert.
    :param fields: Set of fields to add to the FrameData.
    :param frame: Prexisting FrameData to add fields to.
    :return: FrameData with fields requested added from the trajectory.
    """
    if frame is None:
        frame = FrameData()
    if ParticlePositions in fields:
        ParticlePositions.set(
            frame, trajectory.xyz[frame_index] * _MDTrajToNarupa.length
        )
    if ParticleCount in fields:
        frame[ParticleCount] = trajectory.n_atoms
    if BoxVectors in fields and trajectory.unitcell_vectors is not None:
        BoxVectors.set(
            frame, trajectory.unitcell_vectors[frame_index] * _MDTrajToNarupa.length
        )
    mdtraj_topology_to_frame(trajectory.topology, fields=fields, frame=frame)
    return frame


def mdtraj_topology_to_frame(
    topology: Topology,
    *,
    fields: InfiniteSet[str] = MDTRAJ_PROPERTIES,
    frame: Optional[FrameData] = None,
) -> FrameData:
    """
    Convert an MDTraj topology into a FrameData.

    :param topology: Topology to convert.
    :param fields: Set of fields to add to the FrameData.
    :param frame: Prexisting FrameData to add fields to.
    :return: FrameData with fields requested added from the topology.
    """
    if frame is None:
        frame = FrameData()
    if BondPairs in fields:
        BondPairs.set(
            frame, [[bond[0].index, bond[1].index] for bond in topology.bonds]
        )
    if ParticleResidues in fields:
        frame[ParticleResidues] = [atom.residue.index for atom in topology.atoms]
    if ParticleElements in fields:
        frame[ParticleElements] = [atom.element.number for atom in topology.atoms]
    if ParticleNames in fields:
        frame[ParticleNames] = [atom.name for atom in topology.atoms]

    if ResidueNames in fields:
        frame[ResidueNames] = [residue.name for residue in topology.residues]
    if ResidueChains in fields:
        frame[ResidueChains] = [residue.chain.index for residue in topology.residues]

    _get_mdtraj_topology_counts(topology, fields=fields, frame=frame)
    return frame


def _get_mdtraj_topology_counts(
    topology: Topology, *, fields: InfiniteSet[str], frame: FrameData
) -> None:
    if ParticleCount in fields:
        frame[ParticleCount] = topology.n_atoms
    if ResidueCount in fields:
        frame[ResidueCount] = topology.n_residues
    if ChainCount in fields:
        frame[ChainCount] = topology.n_chains
    if BondCount in fields:
        frame[BondCount] = topology.n_bonds
