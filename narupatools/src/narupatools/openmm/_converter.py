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
#
# Originally part of the narupa-openmm package.
# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Modified under the terms of the GPL.

"""Conversion functions between OpenMM objects and Narupa objects."""

import itertools
from typing import Iterable, Optional, Sequence, Type, TypeVar, Union

import numpy as np
from infinite_sets import InfiniteSet, everything
from narupa.trajectory.frame_data import FrameData
from openmm import Context, State, System
from openmm.app import Element, Simulation, Topology

from narupatools.frame import FrameConverter
from narupatools.frame.fields import (
    BondCount,
    BondPairs,
    BoxVectors,
    ChainCount,
    ChainNames,
    KineticEnergy,
    ParticleCount,
    ParticleElements,
    ParticleForces,
    ParticleMasses,
    ParticleNames,
    ParticlePositions,
    ParticleResidues,
    ParticleVelocities,
    PotentialEnergy,
    ResidueChains,
    ResidueCount,
    ResidueIds,
    ResidueNames,
)
from narupatools.openmm._units import UnitsOpenMM
from narupatools.override import override
from narupatools.physics.typing import ScalarArray
from narupatools.physics.units import UnitsNarupa
from narupatools.util import atomic_numbers_to_symbols

DEFAULT_OPENMM_STATE_PROPERTIES = frozenset((ParticlePositions, BoxVectors))

DEFAULT_OPENMM_TOPOLOGY_PROPERTIES = frozenset(
    (
        ResidueNames,
        ResidueChains,
        ResidueCount,
        ParticleCount,
        ChainNames,
        ChainCount,
        ParticleNames,
        ParticleElements,
        ParticleResidues,
        BondPairs,
        BondCount,
        BoxVectors,
    )
)

DEFAULT_OPENMM_SIMULATION_PROPERTIES = frozenset(
    DEFAULT_OPENMM_STATE_PROPERTIES | DEFAULT_OPENMM_TOPOLOGY_PROPERTIES
)

OpenMMToNarupa = UnitsOpenMM >> UnitsNarupa

_TType = TypeVar("_TType")


class OpenMMConverter(FrameConverter):
    """Frame converter for the OpenMM package."""

    @classmethod
    @override(FrameConverter.convert_to_frame)
    def convert_to_frame(  # noqa:D102
        cls,
        object_: _TType,
        /,
        *,
        fields: InfiniteSet[str],
        existing: Optional[FrameData],
    ) -> FrameData:
        if isinstance(object_, Context):
            return openmm_context_to_frame(object_, fields=fields, existing=existing)
        if isinstance(object_, State):
            return openmm_state_to_frame(object_, fields=fields, existing=existing)
        if isinstance(object_, Topology):
            return openmm_topology_to_frame(object_, fields=fields, existing=existing)
        if isinstance(object_, Simulation):
            return openmm_simulation_to_frame(object_, fields=fields, existing=existing)
        raise NotImplementedError

    @classmethod
    @override(FrameConverter.convert_from_frame)
    def convert_from_frame(  # noqa:D102
        cls,
        frame: FrameData,
        destination: Union[Type[_TType], _TType],
        *,
        fields: InfiniteSet[str],
    ) -> _TType:
        if destination == Topology:
            return frame_to_openmm_topology(frame)  # type: ignore
        if isinstance(destination, Simulation):
            copy_frame_to_openmm_simulation(
                frame=frame, simulation=destination, fields=fields
            )
            return destination  # type: ignore
        raise NotImplementedError


def frame_to_openmm_system(frame: FrameData, /) -> System:
    """
    Generate an OpenMM system from a Narupa FrameData.

    :param frame: Narupa FrameData containing masses.
    :raises ValueError: Frame does not have masses to create a system.
    :return: OpenMM system generated using masses.
    """
    system = System()
    if ParticleMasses not in frame:
        raise ValueError("FrameData does not contain masses.")
    for mass in ParticleMasses.get(frame):
        system.addParticle(mass)
    return system


def frame_to_openmm_topology(frame: FrameData, /) -> Topology:
    """
    Convert a Narupa FrameData to an OpenMM topology.

    :param frame: FrameData to convert.
    :return: OpenMM topology populated from the provided FrameData.
    """
    topology = Topology()

    elements = frame[ParticleElements]
    residues = ParticleResidues.get_with_default(frame, itertools.repeat(0))
    resnames = ResidueNames.get_with_default(frame, itertools.repeat("Xxx"))
    names = ParticleNames.get_with_default(frame, atomic_numbers_to_symbols(elements))
    resids = ResidueIds.get_with_default(frame, itertools.repeat(None))
    segs = ResidueChains.get_with_default(frame, itertools.repeat(0))
    segnames = ChainNames.get_with_default(frame, itertools.repeat(None))
    segcount = ChainCount.get_with_default(frame, 1)
    rescount = ResidueCount.get_with_default(frame, 1)
    bonds: Iterable[Sequence[int]] = BondPairs.get_with_default(frame, [])
    box = BoxVectors.get_with_default(frame, None)

    omm_segs = []
    for segname, _ in zip(segnames, range(segcount)):
        omm_segs.append(topology.addChain(id=segname))

    omm_res = []
    for resname, resid, seg, _ in zip(resnames, resids, segs, range(rescount)):
        omm_res.append(topology.addResidue(name=resname, id=resid, chain=omm_segs[seg]))

    omm_atoms = []
    for residue, name, element in zip(residues, names, elements):
        omm_atoms.append(
            topology.addAtom(
                name=name,
                element=Element.getByAtomicNumber(element),
                residue=omm_res[residue],
            )
        )

    for bond in bonds:
        topology.addBond(atom1=omm_atoms[bond[0]], atom2=omm_atoms[bond[1]])

    if box is not None:
        topology.setPeriodicBoxVectors(box)

    return topology


def copy_frame_to_openmm_simulation(
    frame: FrameData, simulation: Simulation, fields: InfiniteSet[str] = everything()
) -> None:
    """Copy fields from a FrameData to a simulation."""
    if ParticlePositions in frame and ParticlePositions in fields:
        simulation.context.setPositions(frame[ParticlePositions])
    if ParticleVelocities in frame and ParticleVelocities in fields:
        simulation.context.setVelocities(frame[ParticleVelocities])


def openmm_simulation_to_frame(
    simulation: Simulation,
    /,
    *,
    fields: InfiniteSet[str] = DEFAULT_OPENMM_SIMULATION_PROPERTIES,
    existing: Optional[FrameData] = None,
) -> FrameData:
    """
    Convert an OpenMM simulation to a Narupa FrameData.

    :param simulation: OpenMM simulation to convert.
    :param fields: Properties to read from frame.
    :param existing: Prexisting FrameData to populate.
    :return: FrameData with requested fields populated from an OpenMM simulation.
    """
    if existing is None:
        frame = FrameData()
    else:
        frame = existing

    if ParticleMasses in fields:
        frame[ParticleMasses] = get_openmm_masses(simulation.system)

    openmm_topology_to_frame(simulation.topology, fields=fields, existing=frame)
    openmm_context_to_frame(simulation.context, fields=fields, existing=frame)

    return frame


def openmm_context_to_frame(
    context: Context,
    /,
    *,
    fields: InfiniteSet[str] = DEFAULT_OPENMM_STATE_PROPERTIES,
    existing: Optional[FrameData] = None,
) -> FrameData:
    """
    Convert an OpenMM context to a Narupa FrameData.

    By converting the context instead of the state, a state can be generated which only
    has the data that is asked for in properties.

    :param context: OpenMM context to convert.
    :param fields: Properties to read from frame.
    :param existing: Prexisting FrameData to populate.
    :return: FrameData with requested fields populated from an OpenMM context.
    """
    if existing is None:
        frame = FrameData()
    else:
        frame = existing

    need_positions = ParticlePositions in fields
    need_forces = ParticleForces in fields
    need_velocities = ParticleVelocities in fields
    need_energy = PotentialEnergy in fields

    state = context.getState(
        getPositions=need_positions,
        getForces=need_forces,
        getEnergy=need_energy,
        getVelocities=need_velocities,
    )

    return openmm_state_to_frame(state, fields=fields, existing=frame)


def openmm_state_to_frame(
    state: State,
    /,
    *,
    fields: InfiniteSet[str] = DEFAULT_OPENMM_STATE_PROPERTIES,
    existing: Optional[FrameData] = None,
) -> FrameData:
    """
    Convert an OpenMM state object to a Narupa FrameData.

    :param state: OpenMM state to convert.
    :param fields: A list of properties to include.
    :param existing: Prexisting FrameData to populate.
    :return: FrameData populated with the requested properties that could be obtained
             from the state object.
    """
    if existing is None:
        frame = FrameData()
    else:
        frame = existing

    if ParticlePositions in fields:
        frame[ParticlePositions] = (
            state.getPositions(asNumpy=True)._value * OpenMMToNarupa.length
        )

    if ParticleForces in fields:
        frame[ParticleForces] = (
            state.getForces(asNumpy=True)._value * OpenMMToNarupa.force
        )

    if ParticleVelocities in fields:
        frame[ParticleVelocities] = (
            state.getVelocities(asNumpy=True)._value * OpenMMToNarupa.velocity
        )

    if PotentialEnergy in fields:
        frame[PotentialEnergy] = (
            state.getPotentialEnergy()._value * OpenMMToNarupa.energy
        )

    if KineticEnergy in fields:
        frame[KineticEnergy] = state.getKineticEnergy()._value * OpenMMToNarupa.energy

    if BoxVectors in fields:
        frame[BoxVectors] = (
            np.array(state.getPeriodicBoxVectors()._value, dtype=float)
            * OpenMMToNarupa.length
        )
    return frame


def openmm_topology_to_frame(
    topology: Topology,
    /,
    *,
    fields: InfiniteSet[str] = DEFAULT_OPENMM_TOPOLOGY_PROPERTIES,
    existing: Optional[FrameData] = None,
) -> FrameData:
    """
    Convert an OpenMM topology object to a Narupa FrameData.

    :param topology: OpenMM topology to convert.
    :param fields: List of properties to include.
    :param existing: Prexisting FrameData to populate.
    :return: FrameData populated with the requested properties that could be obtained
             from the topology object.
    """
    if existing is None:
        frame = FrameData()
    else:
        frame = existing

    if ParticleCount in fields:
        frame[ParticleCount] = topology.getNumAtoms()

    _get_openmm_topology_residue_info(topology, fields=fields, frame=frame)

    if ChainNames in fields:
        frame[ChainNames] = [str(chain.id) for chain in topology.chains()]
    if ChainCount in fields:
        frame[ChainCount] = topology.getNumChains()

    _get_openmm_topology_atom_info(topology, fields=fields, frame=frame)

    _get_openmm_topology_bonds(topology, fields=fields, frame=frame)

    if BoxVectors in fields:
        box = topology.getPeriodicBoxVectors()
        if box is not None:
            frame[BoxVectors] = (
                np.array(box._value, dtype=float) * OpenMMToNarupa.length
            )

    return frame


def _get_openmm_topology_bonds(
    topology: Topology, /, *, fields: InfiniteSet[str], frame: FrameData
) -> None:
    if BondPairs in fields:
        bonds = []
        for bond in topology.bonds():
            bonds.append((bond[0].index, bond[1].index))
        frame[BondPairs] = bonds
    if BondCount in fields:
        frame[BondCount] = topology.getNumBonds()


def _get_openmm_topology_atom_info(
    topology: Topology, /, *, fields: InfiniteSet[str], frame: FrameData
) -> None:
    if (
        ParticleNames in fields
        or ParticleElements in fields
        or ParticleResidues in fields
    ):

        n = topology.getNumAtoms()
        atom_names = np.empty(n, dtype=object)
        elements = np.empty(n, dtype=int)
        residue_indices = np.empty(n, dtype=int)

        for i, atom in enumerate(topology.atoms()):
            atom_names[i] = atom.name
            elements[i] = atom.element.atomic_number
            residue_indices[i] = atom.residue.index

        if ParticleNames in fields:
            frame[ParticleNames] = atom_names
        if ParticleElements in fields:
            frame[ParticleElements] = elements
        if ParticleResidues in fields:
            frame[ParticleResidues] = residue_indices


def _get_openmm_topology_residue_info(
    topology: Topology, /, *, fields: InfiniteSet[str], frame: FrameData
) -> None:
    if ResidueNames in fields:
        frame[ResidueNames] = [residue.name for residue in topology.residues()]
    if ResidueIds in fields:
        frame[ResidueIds] = [residue.id for residue in topology.residues()]
    if ResidueChains in fields:
        frame[ResidueChains] = [residue.chain.index for residue in topology.residues()]
    if ResidueCount in fields:
        frame[ResidueCount] = topology.getNumResidues()


def get_openmm_masses(system: System) -> ScalarArray:
    """
    Get the masses defined in an OpenMM in daltons.

    :param simulation: OpenMM simulation to extract masses from.
    """
    n = system.getNumParticles()
    return np.array([system.getParticleMass(i)._value for i in range(n)], dtype=float)
