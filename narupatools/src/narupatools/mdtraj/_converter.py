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

from narupatools.core.units import UnitsNarupa
from narupatools.frame._converter import FrameConverter
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

MDTRAJ_PROPERTIES = frozenset(
    (
        ParticlePositions.key,
        BondPairs.key,
        ParticleResidues.key,
        ParticleElements.key,
        ParticleNames.key,
        ResidueNames.key,
        ResidueChains.key,
        ParticleCount.key,
        ResidueCount.key,
        BondCount.key,
        ChainCount.key,
        BoxVectors.key,
    )
)

_MDTrajToNarupa = UnitsMDTraj >> UnitsNarupa

_TType = TypeVar("_TType")


class MDTrajConverter(FrameConverter):
    """FrameConverter for the mdtraj package."""

    @classmethod
    @override
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
    @override
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
    if ParticlePositions.key in fields:
        ParticlePositions.set(
            frame, trajectory.xyz[frame_index] * _MDTrajToNarupa.length
        )
    if ParticleCount.key in fields:
        ParticleCount.set(frame, trajectory.n_atoms)
    if BoxVectors.key in fields:
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
    if BondPairs.key in fields:
        BondPairs.set(
            frame, [[bond[0].index, bond[1].index] for bond in topology.bonds]
        )
    if ParticleResidues.key in fields:
        ParticleResidues.set(frame, [atom.residue.index for atom in topology.atoms])
    if ParticleElements.key in fields:
        ParticleElements.set(frame, [atom.element.number for atom in topology.atoms])
    if ParticleNames.key in fields:
        ParticleNames.set(frame, [atom.name for atom in topology.atoms])

    if ResidueNames.key in fields:
        ResidueNames.set(frame, [residue.name for residue in topology.residues])
    if ResidueChains.key in fields:
        ResidueChains.set(frame, [residue.chain.index for residue in topology.residues])

    _get_mdtraj_topology_counts(topology, fields=fields, frame=frame)
    return frame


def _get_mdtraj_topology_counts(
    topology: Topology, *, fields: InfiniteSet[str], frame: FrameData
) -> None:
    if ParticleCount.key in fields:
        ParticleCount.set(frame, topology.n_atoms)
    if ResidueCount.key in fields:
        ResidueCount.set(frame, topology.n_residues)
    if ChainCount.key in fields:
        ChainCount.set(frame, topology.n_chains)
    if BondCount.key in fields:
        BondCount.set(frame, topology.n_bonds)
