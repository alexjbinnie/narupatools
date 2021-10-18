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
# Originally part of the narupa-ase package.
# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Modified under the terms of the GPL.

"""Conversion functions for converting ASE objects to Narupa objects."""

import contextlib
import itertools
from contextlib import suppress
from typing import Any, Dict, List, Optional, Type, TypeVar, Union

import numpy as np
from ase.atoms import Atoms
from ase.calculators.calculator import PropertyNotImplementedError
from infinite_sets import InfiniteSet, everything
from narupa.trajectory.frame_data import FrameData

from narupatools.ase._units import UnitsASE
from narupatools.ase.calculators import ConstantCalculator
from narupatools.frame import FrameConverter
from narupatools.frame.fields import (
    BondPairs,
    BondTypes,
    BoxVectors,
    KineticEnergy,
    ParticleCharges,
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
    ResidueNames, ChainCount,
)
from narupatools.frame.util import calculate_residue_entities
from narupatools.mdanalysis import UnitsMDAnalysis
from narupatools.override import override
from narupatools.physics.units import UnitsNarupa

_ASEToNarupa = UnitsASE >> UnitsNarupa
_NarupaToASE = UnitsNarupa >> UnitsASE

_MDAnalysisToASE = UnitsMDAnalysis >> UnitsASE

ASE_PROPERTIES = frozenset(
    (
        ParticlePositions,
        ParticleElements,
        BoxVectors,
        ParticleCount,
        ResidueNames,
        ParticleNames,
        ParticleResidues,
        ResidueCount,
        BondPairs,
        BondTypes,
        ResidueChains,
    )
)

_TType = TypeVar("_TType")


class ASEConverter(FrameConverter):
    """Converters for the ASE package."""

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
        if isinstance(object_, Atoms):
            return ase_atoms_to_frame(object_, fields=fields)
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
        if destination == Atoms:
            return frame_to_ase_atoms(frame=frame, fields=fields)  # type: ignore
        if isinstance(destination, Atoms):
            copy_frame_to_ase_atoms(atoms=destination, frame=frame, fields=fields)
            return destination
        raise NotImplementedError


def copy_frame_to_ase_atoms(
    *, atoms: Atoms, frame: FrameData, fields: InfiniteSet[str] = everything()
) -> None:
    if ParticlePositions in fields and ParticlePositions in frame:
        atoms.set_positions(frame[ParticlePositions] * _NarupaToASE.length)
    if ParticleVelocities in fields and ParticleVelocities in frame:
        atoms.set_velocities(frame[ParticleVelocities] * _NarupaToASE.velocity)
    if ParticleMasses in fields and ParticleMasses in frame:
        atoms.set_masses(frame[ParticleMasses] * _NarupaToASE.masses)


def frame_to_ase_atoms(
    frame: FrameData, *, fields: InfiniteSet[str] = everything()
) -> Atoms:
    """
    Convert a Narupa FrameData to an ASE atoms object.

    :param frame: Narupa FrameData to convert.
    :param fields: Set of fields to convert.
    :return: ASE atoms objects with properties read from the frame.
    """
    kwargs: Dict[str, Any] = {}
    if ParticlePositions in fields:
        with contextlib.suppress(KeyError):
            kwargs["positions"] = ParticlePositions.get(frame) * _NarupaToASE.length
    if ParticleMasses in fields:
        with contextlib.suppress(KeyError):
            kwargs["masses"] = ParticleMasses.get(frame) * _NarupaToASE.mass
    if ParticleVelocities in fields:
        with contextlib.suppress(KeyError):
            kwargs["velocities"] = ParticleVelocities.get(frame) * _NarupaToASE.velocity
    if ParticleElements in fields:
        with contextlib.suppress(KeyError):
            kwargs["numbers"] = ParticleElements.get(frame)
    if BoxVectors in fields:
        with contextlib.suppress(KeyError):
            kwargs["cell"] = BoxVectors.get(frame) * _NarupaToASE.length
    if ParticleCharges in fields:
        with contextlib.suppress(KeyError):
            kwargs["charges"] = ParticleCharges.get(frame) * _NarupaToASE.charge

    calc_kwargs: Dict[str, Any] = {}

    if ParticleForces in fields:
        with contextlib.suppress(KeyError):
            calc_kwargs["forces"] = ParticleForces.get(frame) * _NarupaToASE.force
    if PotentialEnergy in fields:
        with contextlib.suppress(KeyError):
            calc_kwargs["energy"] = PotentialEnergy.get(frame) * _NarupaToASE.energy
    if ParticleCharges in fields:
        with contextlib.suppress(KeyError):
            calc_kwargs["charges"] = ParticleCharges.get(frame) * _NarupaToASE.charge

    calculator = ConstantCalculator(**calc_kwargs)

    atoms = Atoms(**kwargs)

    _add_bonds_to_ase_atoms(frame, fields, atoms)

    if ParticleResidues in fields and ParticleResidues in frame:
        particle_residues = ParticleResidues.get(frame)
        with contextlib.suppress(KeyError):
            atoms.set_array("residuenumbers", ParticleResidues.get(frame))

        if ResidueNames in fields and ResidueNames in frame:
            residue_names = ResidueNames.get(frame)
            atoms.set_array(
                "residuenames", np.array([residue_names[i] for i in particle_residues])
            )

    if ParticleNames in fields and ParticleNames in frame:
        atoms.set_array("atomtypes", np.array(ParticleNames.get(frame)))

    atoms.set_calculator(calculator)
    return atoms


def ase_atoms_to_frame(
    atoms: Atoms,
    *,
    fields: InfiniteSet[str] = ASE_PROPERTIES,
    frame: Optional[FrameData] = None,
) -> FrameData:
    """
    Convert an ASE Atoms object to a Narupa FrameData.

    :param atoms: ASE Atoms object to convert.
    :param fields: A collection ofs that should be added to the frame if available.
    :param frame: An optional preexisting FrameData to populate.
    :return: A FrameData populated with information available in the ASE atoms object
             whoses are present in the properties parameter.
    """
    if frame is None:
        frame = FrameData()

    _add_ase_atoms_particles_to_frame(atoms, fields, frame)
    _add_ase_atoms_residues_to_frame(atoms, fields, frame)
    _add_ase_atoms_calculated_properties_to_frame(atoms, fields, frame)

    if BoxVectors in fields:
        frame[BoxVectors] = atoms.get_cell() * _ASEToNarupa.length

    if KineticEnergy in fields:
        frame[KineticEnergy] = atoms.get_kinetic_energy() * _ASEToNarupa.energy

    _add_ase_bonds_to_frame(atoms, fields, frame)

    if (
        BondPairs in frame
        and ParticleResidues in frame
        and ResidueCount in frame
        and ResidueChains not in frame
        and len(frame[BondPairs]) > 0
        and ResidueChains in fields
    ):
        frame[ResidueChains] = calculate_residue_entities(
            residue_count=frame[ResidueCount],
            particle_residues=frame[ParticleResidues],
            bond_pairs=frame[BondPairs],
        )
        frame[ChainCount] = frame[ResidueChains].max() + 1

    return frame


def _add_ase_atoms_particles_to_frame(
    atoms: Atoms, fields: InfiniteSet[str], frame: FrameData
) -> None:
    if ParticlePositions in fields:
        frame[ParticlePositions] = atoms.get_positions() * _ASEToNarupa.length

    if ParticleCount in fields:
        frame[ParticleCount] = len(atoms)

    if ParticleElements in fields:
        elements = []
        for atom in atoms:
            elements.append(atom.number)

        frame[ParticleElements] = elements

    if ParticleMasses in fields:
        frame[ParticleMasses] = atoms.get_masses()

    if ParticleVelocities in fields:
        frame[ParticleVelocities] = atoms.get_velocities() * _ASEToNarupa.velocity

    if ParticleNames in fields:
        if "atomtypes" in atoms.arrays:
            frame[ParticleNames] = atoms.arrays["atomtypes"]
        else:
            frame[ParticleNames] = list(atoms.symbols)


def _add_ase_atoms_calculated_properties_to_frame(
    atoms: Atoms, fields: InfiniteSet[str], frame: FrameData
) -> None:
    if ParticleCharges in fields:
        calc_charges = False
        if atoms.calc is not None:
            with contextlib.suppress(PropertyNotImplementedError):
                frame[ParticleCharges] = atoms.get_charges()
                calc_charges = True
        if not calc_charges:
            frame[ParticleCharges] = atoms.get_initial_charges()

    if ParticleForces in fields and atoms.calc is not None:
        with contextlib.suppress(PropertyNotImplementedError):
            frame[ParticleForces] = atoms.get_forces() * _ASEToNarupa.force

    if PotentialEnergy in fields and atoms.calc is not None:
        with suppress(PropertyNotImplementedError):
            PotentialEnergy.set(
                frame, atoms.get_potential_energy() * _ASEToNarupa.energy
            )


def _add_ase_atoms_residues_to_frame(
    atoms: Atoms, fields: InfiniteSet[str], frame: FrameData
) -> None:
    if (
        ResidueNames in fields or ResidueCount in fields or ParticleResidues in fields
    ) and "residuenumbers" in atoms.arrays:
        segid_to_index: Dict[Any, int] = {}
        res_to_first_particle_index = []
        index = 0
        for atom_index, segid in enumerate(atoms.arrays["residuenumbers"]):
            if segid not in segid_to_index:
                segid_to_index[segid] = index
                res_to_first_particle_index.append(atom_index)
                index += 1

        if ParticleResidues in fields:
            ParticleResidues.set(
                frame,
                [segid_to_index[segid] for segid in atoms.arrays["residuenumbers"]],
            )

        if ResidueNames in fields and "residuenames" in atoms.arrays:
            ResidueNames.set(
                frame,
                [
                    str(atoms.arrays["residuenames"][atom_index]).strip()
                    for atom_index in res_to_first_particle_index
                ],
            )

        if ResidueCount in fields and "residuenames" in atoms.arrays:
            frame[ResidueCount] = len(res_to_first_particle_index)


def get_bonds(atoms: Atoms):
    bonds = []
    for index, atom_bonds in enumerate(atoms.get_array("bonds", copy=False)):
        for other in atom_bonds:
            bonds.append([index, other])
    return np.array(bonds)


def set_bonds(atoms: Atoms, src_bonds: Any, /):
    bonds: List[List[int]] = [[] for _ in range(len(atoms))]

    for bond in src_bonds:
        i = min(bond[0], bond[1])
        bonds[i].append(max(bond[0], bond[1]))

    atoms.arrays["bonds"] = np.array(bonds, dtype=object)


Atoms.get_bonds = get_bonds
Atoms.set_bonds = set_bonds


def _add_bonds_to_ase_atoms(
    frame: FrameData, fields: InfiniteSet[str], atoms: Atoms
) -> None:

    bonds: List[List[int]] = [[] for _ in range(len(atoms))]
    bond_types: List[List[str]] = [[] for _ in range(len(atoms))]

    src_bonds: Any = BondPairs.get_with_default(frame, [])
    src_types = BondTypes.get_with_default(frame, itertools.repeat("s"))

    for bond, bond_type in zip(src_bonds, src_types):
        i = min(bond[0], bond[1])
        bonds[i].append(max(bond[0], bond[1]))
        bond_types[i].append(bond_type)

    if BondPairs in fields:
        atoms.arrays["bonds"] = np.array(bonds, dtype=object)
    if BondTypes in fields:
        atoms.arrays["bond_types"] = np.array(bond_types, dtype=object)


def _add_ase_bonds_to_frame(
    atoms: Atoms, fields: InfiniteSet[str], frame: FrameData
) -> None:

    if BondPairs in fields and atoms.has("bonds"):
        bonds = []
        for index, atom_bonds in enumerate(atoms.get_array("bonds", copy=False)):
            for other in atom_bonds:
                bonds.append([index, other])
        frame[BondPairs] = bonds
    if BondTypes in fields and atoms.has("bond_types"):
        bond_types = []
        for types in atoms.get_array("bond_types", copy=False):
            for type_ in types:
                bond_types.append(type_)
        frame[BondTypes] = bond_types
