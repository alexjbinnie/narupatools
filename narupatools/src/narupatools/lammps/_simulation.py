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

"""Wrapper around LAMMPS simulation."""

from __future__ import annotations

import contextlib
import ctypes
from abc import ABCMeta, abstractmethod
from threading import Lock
from typing import Any, Dict, Literal, Optional, Set, TypeVar, Union, overload

import numpy as np
import numpy.typing as npt
import quaternion
from infinite_sets import InfiniteSet, everything
from narupa.trajectory import FrameData

import narupatools.lammps.atom_properties as PROPERTIES
import narupatools.lammps.computes as COMPUTES
import narupatools.lammps.globals as GLOBALS
import narupatools.lammps.settings as SETTINGS
from narupatools.frame import FrameSource
from narupatools.frame.fields import (
    BondPairs,
    BoxPeriodic,
    BoxVectors,
    KineticEnergy,
    ParticleCharges,
    ParticleCount,
    ParticleElements,
    ParticleForces,
    ParticleMasses,
    ParticlePositions,
    ParticleResidues,
    ParticleVelocities,
    PotentialEnergy,
    ResidueCount,
)
from narupatools.lammps._units import get_unit_system
from narupatools.physics.energy import kinetic_energy
from narupatools.physics.transformation import Rotation
from narupatools.physics.typing import Vector3
from narupatools.physics.units import UnitsNarupa, UnitSystem, degree, radian
from narupatools.physics.vector import magnitude, normalized, vector
from narupatools.util import mass_to_element

from ._constants import VariableDimension, VariableType
from ._exception_wrapper import catch_lammps_warnings_and_exceptions
from ._wrapper import Extractable, LAMMPSWrapper
from .exceptions import UnknownAtomPropertyError
from .regions import Region, RegionSpecification

_TReturnType = TypeVar("_TReturnType")


class LAMMPSSimulation(FrameSource):
    """Wrapper around a LAMMPS simulation."""

    def __init__(self, lammps: LAMMPSWrapper, units: Optional[UnitSystem] = None):
        self.__lammps = lammps

        unit_str = self.extract(GLOBALS.UnitStyle)
        if unit_str == "lj":
            if units is None:
                raise ValueError(
                    "LAMMPS simulation uses 'lj' units but no UnitSystem provided."
                )
            self.unit_system = units
        else:
            if units is not None:
                raise ValueError(
                    "LAMMPS simulation does not use 'lj' units but UnitSystem provided."
                )
            self.unit_system = get_unit_system(unit_str)
        self._lammps_to_narupa = self.unit_system >> UnitsNarupa
        self._narupa_to_lammps = UnitsNarupa >> self.unit_system
        self.indexing = LAMMPSIndexing(self)
        self.bond_compute: Optional[COMPUTES.ComputeReference] = None

        self._run_lock = Lock()

        self._needs_pre_run = True

        self._mass_compute = COMPUTES.AtomProperty(
            self.__lammps,
            compute_id="mass",
            properties=["mass"],
            datatype=VariableType.DOUBLE,
        )

        self._kinetic_energy_compute = COMPUTES.KineticEnergy.create(
            self.__lammps, compute_id="ke"
        )

        self._potential_energy_compute = COMPUTES.ComputeGlobalReference.create(
            self.__lammps, compute_id="thermo_pe", dimension=VariableDimension.SCALAR
        )

        self._temperature_compute = COMPUTES.ComputeGlobalReference.create(
            self.__lammps, compute_id="thermo_temp", dimension=VariableDimension.SCALAR
        )
        self._pressure_compute = COMPUTES.ComputeGlobalReference.create(
            self.__lammps, compute_id="thermo_press", dimension=VariableDimension.SCALAR
        )

        self.atoms = AtomArray(self)
        self.types = TypeArray(self)
        self.mols = MolArray(self)

        self._quat_compute: Optional[
            COMPUTES.ComputeAtomReference[npt.NDArray[np.float64]]
        ] = None
        self._angmom_compute: Optional[
            COMPUTES.ComputeAtomReference[npt.NDArray[np.float64]]
        ] = None

    @property
    def _wrapper(self) -> LAMMPSWrapper:
        return self.__lammps

    @classmethod
    def create_new(
        cls,
        units: Union[
            UnitSystem,
            Literal["real", "metal", "si", "cgs", "electron", "micro", "nano"],
        ],
    ) -> LAMMPSSimulation:
        """
        Create a new LAMMPS simulation with the given units system.

        This automatically runs the 'atom_modify map yes' command required to use atom
        IDs and hence allow certain operations such as setting positions.
        :param units: LAMMPS unit system to use.
        :return: LAMMPS simulation in provided units.
        """
        lammps = LAMMPSWrapper()
        with catch_lammps_warnings_and_exceptions():
            lammps.command("atom_modify map yes")
            if isinstance(units, UnitSystem):
                lammps.command("units lj")
            else:
                lammps.command(f"units {units}")
        if isinstance(units, UnitSystem):
            return cls(lammps, units)
        else:
            return cls(lammps)

    @classmethod
    def from_file(
        cls, filename: str, units: Optional[UnitSystem] = None
    ) -> LAMMPSSimulation:
        """
        Load LAMMPS simulation from a file.

        This automatically runs the 'atom_modify map yes' command required to use atom
        IDs and hence allow certain operations such as setting positions.
        :param filename: Filename of LAMMPS input.
        :param units: If the file uses LJ units, a UnitSystem must be provided.
        :return: LAMMPS simulation based on file.
        """
        lammps = LAMMPSWrapper()
        with catch_lammps_warnings_and_exceptions():
            lammps.command("atom_modify map yes")
            lammps.file(filename)
            lammps.command("run 0")
        return cls(lammps, units)

    def __len__(self) -> int:
        return len(self.__lammps)

    def gather_atoms(
        self, atom_property: PROPERTIES.AtomProperty[_TReturnType]
    ) -> _TReturnType:
        """
        Gather an atom property from across all processors, and order by Atom ID.

        :param atom_property: The atom property to gather.
        :return: The atom data as a numpy array.
        """
        return self.__lammps.gather_atoms(  # type: ignore[return-value]
            atom_property.key, atom_property.components
        )

    def scatter_atoms(
        self, atom_property: PROPERTIES.AtomProperty, value: np.ndarray
    ) -> None:
        """
        Scatter an atom property across all processors.

        :param atom_property: The atom property to scatter.
        :param value: The value to distribute.
        """
        self.__lammps.scatter_atoms(
            atom_property.key, atom_property.datatype, atom_property.components, value
        )

    def run(self, steps: int = 1, pre: bool = False) -> None:
        """
        Run the simulation for a given number of steps.

        :param steps: Number of steps to run.
        """
        with self._run_lock:
            should_pre = "yes" if (self._needs_pre_run or pre) else "no"
            self.__lammps.command(f"run {steps} pre {should_pre} post no")
            self._needs_pre_run = False

    @property
    def potential_energy(self) -> float:
        """
        Potential energy of the system in current LAMMPS units.

        The potential energy is computed using the 'thermo_pe' compute which is added to
        all LAMMPS simulations automatically. This uses the 'pe' compute style which
        includes contributions from pairs, bonds, angles, dihedrals, impropers, kspace
        and various fixes.
        """
        return self._potential_energy_compute.extract()

    @property
    def kinetic_energy(self) -> float:
        """Kinetic energy of the system in current LAMMPS units."""
        return self._kinetic_energy_compute.extract()

    @property
    def temperature(self) -> float:
        r"""
        Temperature of the system in current LAMMPS units.

        The temperature is computed using the 'thermo_temp' compute which is added to
        all LAMMPS simulations automatically. This uses the 'temp' compute style which
        uses the equation:

        .. math:: KE = \frac{D}{2} N k T

        where :math:`KE` is the kinetic energy, :math:`D` is the dimension of the
        system, :math:`N` is the number of atoms and :math:`k` is the Boltzmann
        constant.

        This removes degrees of freedom due to other fixes that limit molecular motion.
        """
        return self._temperature_compute.extract()

    @property
    def pressure(self) -> float:
        """Pressure of the system in current LAMMPS units."""
        return self._pressure_compute.extract()

    @property
    def forces(self) -> npt.NDArray[np.float64]:
        """Forces on each atom in current LAMMPS units."""
        return self.gather_atoms(PROPERTIES.Force)

    @property
    def atom_ids(self) -> npt.NDArray[np.int64]:
        """Atom ID of each atom as defined by LAMMPS."""
        return self.gather_atoms(PROPERTIES.AtomID)

    @property
    def charges(self) -> npt.NDArray[np.float64]:
        """Charge of each atom in current LAMMPS units."""
        try:
            return self.gather_atoms(PROPERTIES.Charge)
        except UnknownAtomPropertyError as e:
            raise AttributeError from e

    @property
    def masses(self) -> npt.NDArray[np.float64]:
        """Masses of each atom in current LAMMPS units."""
        return self._mass_compute.gather()

    @property
    def positions(self) -> npt.NDArray[np.float64]:
        """Positions of each atom in current LAMMPS units."""
        return self.gather_atoms(PROPERTIES.Position)

    @positions.setter
    def positions(self, positions: np.ndarray) -> None:
        self.scatter_atoms(PROPERTIES.Position, positions)

    @property
    def velocities(self) -> npt.NDArray[np.float64]:
        """Velocities of each atom in current LAMMPS units."""
        return self.gather_atoms(PROPERTIES.Velocity)

    @velocities.setter
    def velocities(self, velocities: np.ndarray) -> None:
        self.scatter_atoms(PROPERTIES.Velocity, velocities)

    @property
    def orientations(self) -> npt.NDArray[quaternion.quaternion]:
        """Orientations of each atom as unit quaternions."""
        if self._quat_compute is None:
            if not self.extract(SETTINGS.AtomStylesCanBeEllipsoid):
                raise AttributeError("orientations not supported for this simulation.")
            self._quat_compute = COMPUTES.AtomProperty(
                self.__lammps,
                compute_id="quat",
                properties=["quatw", "quati", "quatj", "quatk"],
                datatype=VariableType.DOUBLE,
            )
        return quaternion.as_quat_array(self._quat_compute.gather())

    @property
    def ellipsoid_axes(self) -> npt.NDArray[np.float64]:
        """Principle moments of inertia of each atom ."""
        if not hasattr(self, "_shape_compute"):
            if not self.extract(SETTINGS.AtomStylesCanBeEllipsoid):
                raise AttributeError(
                    "ellipsoid_axes not supported for this simulation."
                )
            self._shape_compute = COMPUTES.AtomProperty(
                self.__lammps,
                compute_id="shape",
                properties=["shapex", "shapey", "shapez"],
                datatype=VariableType.DOUBLE,
            )
        return self._shape_compute.gather()

    @property
    def angular_momenta(self) -> npt.NDArray[np.float64]:
        """Angular momentum of each particle."""
        if self._angmom_compute is None:
            if not self.extract(SETTINGS.AtomStylesCanBeEllipsoid):
                raise AttributeError(
                    "angular momenta not supported for this simulation."
                )
            self._angmom_compute = COMPUTES.AtomProperty(
                self.__lammps,
                compute_id="angmom",
                properties=["angmomx", "angmomy", "angmomz"],
                datatype=VariableType.DOUBLE,
            )
        return self._angmom_compute.gather()

    @property
    def bond_energies(self) -> npt.NDArray[np.float64]:
        """Get the energy of each bond."""
        if not hasattr(self, "_bond_energy_compute"):
            self._bond_energy_compute = COMPUTES.BondLocal(
                self.__lammps, compute_id="bond_energy", properties=["engpot"]
            )
        return self._bond_energy_compute.extract()

    def setup_imd(self, callback: Any) -> None:
        """
        Setup IMD for this simulation.

        The callback must take no parameters and return a tuple of the energy, forces and torques, all in
        internal LAMMPS units. Both the forces and the torques may be zero if they do not exist.
        """
        self._imd_calculate = callback
        self._setup_imd_forces()

    def _setup_imd_forces(self) -> None:
        self.command("fix imd_force all external pf/callback 1 1")
        self._wrapper._lammps.set_fix_external_callback(  # type: ignore
            "imd_force", self._imd_callback
        )

    def _imd_callback(
        self,
        lmp: Any,
        ntimestep: int,
        nlocal: int,
        ids: np.ndarray,
        x: np.ndarray,
        f: np.ndarray,
    ) -> None:
        """Invoked by the IMD fix external to apply forces."""
        if lmp is None:
            lmp = self._wrapper._lammps
        energy, forces, torques = self._imd_calculate()
        f[:, :] = 0
        if forces is not None:
            f[:] = forces[ids - 1]
        lmp.fix_external_set_energy_global("imd_force", energy)
        if torques is not None:
            lmp.scatter_atoms_subset(
                "torque",
                VariableType.DOUBLE_ARRAY,
                3,
                len(ids),
                ids.ctypes.data_as(ctypes.POINTER(lmp.c_tagint)),
                torques[ids - 1].ctypes.data_as(ctypes.c_void_p),
            )

    def command(self, command: str) -> None:
        """
        Run an arbitrary LAMMPS command.

        :param command: LAMMPS command to run.
        """
        self.__lammps.command(command)
        self._needs_pre_run = True

    def get_frame(
        self,
        *,
        fields: InfiniteSet[str] = everything(),
        existing: Optional[FrameData] = None,
    ) -> FrameData:
        """
        Create a Narupa FrameData of the LAMMPS simulation.

        :param fields: Fields to populate of the frame data.
        :param existing: Preexisting frame data.
        :return: FrameData with the appropriate fields.
        """
        if existing is None:
            frame: FrameData = FrameData()
        else:
            frame = existing

        if self._run_lock:
            if ParticleCount in fields:
                frame[ParticleCount] = len(self)
            if ParticlePositions in fields:
                frame[ParticlePositions] = (
                    self.positions * self._lammps_to_narupa.length
                )
            if ParticleVelocities in fields:
                frame[ParticleVelocities] = (
                    self.velocities * self._lammps_to_narupa.velocity
                )
            if ParticleMasses in fields:
                frame[ParticleMasses] = self.masses * self._lammps_to_narupa.mass
            if ParticleElements in fields:
                elements = np.vectorize(mass_to_element)(
                    self.masses * self._lammps_to_narupa.mass
                )
                if not np.any(np.equal(elements, None)):  # type: ignore[call-overload]
                    frame[ParticleElements] = elements
            if ParticleCharges in fields:
                with contextlib.suppress(AttributeError):
                    frame[ParticleCharges] = (
                        self.charges * self._lammps_to_narupa.charge
                    )
            if ParticleForces in fields:
                frame[ParticleForces] = self.forces * self._lammps_to_narupa.force
            if PotentialEnergy in fields:
                frame[PotentialEnergy] = (
                    self.potential_energy * self._lammps_to_narupa.energy
                )
            if KineticEnergy in fields:
                frame[KineticEnergy] = kinetic_energy(
                    velocities=self.velocities, masses=self.masses
                )
            if fields & {BoxVectors, BoxPeriodic}:
                box, periodic = self.__lammps.extract_box_vectors()
                if BoxVectors in fields:
                    frame[BoxVectors] = box * self._lammps_to_narupa.length
                if BoxPeriodic in fields:
                    frame[BoxPeriodic] = periodic

            if (
                BondPairs in fields
                and self[SETTINGS.AtomStylesIncludesMolecularTopology]
            ):
                self.indexing.recompute()
                if self.bond_compute is None:
                    self.bond_compute = COMPUTES.LocalProperty(
                        self.__lammps,
                        compute_id="bonds",
                        properties=["btype", "batom1", "batom2"],
                    )
                bonds = self.bond_compute.extract()
                frame[BondPairs] = self.indexing.atom_id_to_ordered(bonds[:, 1:])

            if (ParticleResidues in fields or ResidueCount in fields) and self[
                SETTINGS.AtomStylesIncludesMolecularTopology
            ]:
                self.indexing.recompute()
                mol_ids = self.gather_atoms(PROPERTIES.MoleculeID)

                unique_mol_ids = np.unique(mol_ids)
                mol_id_to_index = {
                    mol_id: mol_index
                    for mol_index, mol_id in zip(
                        np.argsort(unique_mol_ids), unique_mol_ids
                    )
                }
                frame[ParticleResidues] = np.vectorize(mol_id_to_index.get)(mol_ids)
                frame[ResidueCount] = len(unique_mol_ids)

        return frame

    def __getitem__(self, item: Extractable[_TReturnType]) -> _TReturnType:
        return item.extract(self.__lammps)

    @property
    def timestep(self) -> float:
        """Timestep of the simulation in current LAMMPS units."""
        return self.extract(GLOBALS.TimestepLength)

    @timestep.setter
    def timestep(self, value: float) -> None:
        self.command(f"timestep {value}")
        self._timestep = None

    def create_atom(
        self, *, atom_type: int, position: Vector3, rotation: Optional[Rotation] = None
    ) -> SingleAtomReference:
        """
        Insert an atom into the simulation.

        :param atom_type: Type of the atom to insert.
        :param position: Position of the atom in current LAMMPS units.
        :param rotation: Initial rotation of the atom.
        """
        if len(self) == 0:
            prev_ids: Set[np.int64] = set()
        else:
            prev_ids = set(self.gather_atoms(PROPERTIES.AtomID))  # type: ignore[arg-type]
        position = position
        command = (
            f"create_atoms {atom_type} single {position[0]} {position[1]} {position[2]}"
        )
        if rotation is not None:
            rot_vec = rotation.rotation_vector
            theta = magnitude(rot_vec) * (radian >> degree)
            axis = normalized(rot_vec)
            if theta == 0:
                axis = vector(1, 0, 0)
            command += f" rotate {theta} {axis[0]} {axis[1]} {axis[2]}"
        self.command(command)
        new_ids: Set[np.int64] = set(self.gather_atoms(PROPERTIES.AtomID))  # type: ignore[arg-type]
        diff = new_ids - prev_ids
        assert len(diff) == 1
        return self.atoms[int(list(diff)[0])]

    def create_region(
        self, region: RegionSpecification, *, region_id: Optional[str] = None
    ) -> Region:
        """
        Create a region in the simulation.

        :param region: Specification to the region.
        :param region_id: Optional id for the region.
        :return: Region object.
        """
        if region_id is None:
            region_id = self.__lammps.regions.generate_id()
        self.command(f"region {region_id} {region.style} {region.args}")
        return Region(self, region_id, region)

    def create_box(
        self, n_types: int, region: Union[Region, RegionSpecification]
    ) -> None:
        """
        Create an empty simulation box with the given number of atom types.

        :param n_types: Number of atom types.
        :param region: Region to define.
        :raises ValueError: Region does belong to this simulation.
        """
        if isinstance(region, RegionSpecification):
            region = self.create_region(region)
        elif region._simulation is not self:
            raise ValueError("Region does not belong to this simulation!")
        self.command(f"create_box {n_types} {region.region_id}")

    def setup_langevin(self, *, temperature: float, friction: float, seed: int) -> None:
        """
        Setup a Langevin thermostat.

        :param temperature: Temperature in LAMMPS units
        :param friction: Friction of the thermostat in LAMMPS units.
        :param seed: Seed to use for random number generation.
        """
        self.command(
            f"fix langevin all langevin {temperature} {temperature} {1.0 / friction} "
            f"{seed}"
        )

    def __del__(self) -> None:
        self.__lammps.close()

    def file(self, filename: str) -> None:
        """Runn all commands found in the provided input file."""
        with catch_lammps_warnings_and_exceptions():
            self.__lammps.file(filename)

    def extract(self, obj: Extractable[_TReturnType]) -> _TReturnType:
        """Extract a value from the simulation."""
        return obj.extract(self.__lammps)


class LAMMPSIndexing:
    """
    Manages mappings between different ways that atoms are indexed.

    There are three important orderings to consider when using LAMMPS with narupatools. They are:

    * The actual ordering of the atoms that LAMMPS uses internally.
    * Each atom has a unique integer Atom ID.
    * The ordering that narupatools exposes the atoms as.

    Importantly, using gather_atoms orders the atoms by Atom ID, whilst extract_compute does not.
    """

    def __init__(self, simulation: LAMMPSSimulation):
        self._simulation = simulation
        self.__atom_ids = COMPUTES.AtomProperty(
            simulation._wrapper,
            compute_id="atomids",
            properties=["id"],
            datatype=VariableType.INTEGER,
        )
        self._unordered_to_atom_ids: np.ndarray = np.array([])
        self._ordered_to_unordered: np.ndarray = np.array([])
        self._ordered_to_atom_ids: np.ndarray = np.array([])
        self._atom_ids_to_ordered: Dict[int, int] = {}

    def recompute(self) -> None:
        """Recompute the indexing of atom ids."""
        self._unordered_to_atom_ids = self.__atom_ids.extract()
        self._ordered_to_unordered = np.argsort(self._unordered_to_atom_ids)
        self._ordered_to_atom_ids = self._unordered_to_atom_ids[
            self._ordered_to_unordered
        ]
        for ordered_index, atom_id in enumerate(self._ordered_to_atom_ids):
            self._atom_ids_to_ordered[atom_id] = ordered_index

    def unordered_to_ordered(self, array: np.ndarray) -> np.ndarray:
        """Take an array ordered by LAMMPS and sort it by Atom ID."""
        return array[self._ordered_to_unordered]  # type: ignore[no-any-return]

    def atom_id_to_ordered(self, array: np.ndarray) -> np.ndarray:
        """Maps Atom IDs to zero-based ordered index."""
        return np.vectorize(self._atom_ids_to_ordered.get)(array)  # type: ignore[no-any-return]

    def ordered_to_atom_id(self, index: int) -> int:
        """Maps a zero-based index to an Atom ID."""
        return self._ordered_to_atom_ids[index]  # type: ignore[no-any-return]


def _slice_to_lammps(slice_: slice) -> str:
    """Convert a python slice to a LAMMPS representation."""
    if slice_.start is None:
        if slice_.stop is None:
            return "*"
        else:
            return f"*{slice_.stop}"
    else:
        if slice_.stop is None:
            return f"{slice_.start}*"
        else:
            return f"{slice_.start}*{slice_.stop}"


class AtomArray:
    """Structure allowing access to atoms of a simulation by ID."""

    def __init__(self, simulation: LAMMPSSimulation):
        self._simulation = simulation

    @overload
    def __getitem__(self, index: int, /) -> SingleAtomReference:
        ...

    @overload
    def __getitem__(self, index: slice, /) -> AtomRangeReference:
        ...

    def __getitem__(self, index: Union[int, slice]) -> AtomSetReference:
        if isinstance(index, (int, np.integer)):
            return SingleAtomReference(self._simulation, int(index))
        else:
            return AtomRangeReference(self._simulation, index)


class TypeArray:
    """Structure allowing access to atom types of a simulation by ID."""

    def __init__(self, simulation: LAMMPSSimulation):
        self._simulation = simulation

    def __getitem__(
        self, index: Union[int, slice]
    ) -> Union[SingleTypeReference, TypeRangeReference]:
        if isinstance(index, int):
            return SingleTypeReference(self._simulation, index)
        else:
            return TypeRangeReference(self._simulation, index)


class MolArray:
    """Structure allowing access to mols of a simulation by ID."""

    def __init__(self, simulation: LAMMPSSimulation):
        self._simulation = simulation

    def __getitem__(
        self, index: Union[int, slice]
    ) -> Union[SingleMolReference, MolRangeReference]:
        if isinstance(index, int):
            return SingleMolReference(self._simulation, index)
        else:
            return MolRangeReference(self._simulation, index)


class AtomSetReference(metaclass=ABCMeta):
    """Base class for a reference to a set of atoms."""

    def __init__(self, simulation: LAMMPSSimulation):
        self._simulation = simulation

    @abstractmethod
    def _selection_str(self) -> str:
        ...

    @abstractmethod
    def _type_str(self) -> str:
        ...

    def _set(self, command: str) -> None:
        self._simulation.command(
            f"set {self._type_str()} {self._selection_str()} {command}"
        )

    def set_type(self, *, atom_type: int) -> None:
        """
        Set the atom type of the given atom(s).

        :param atom_type: Atom type.
        """
        self._set(f"type {atom_type}")

    def set_quat(self, *, rotation: Rotation) -> None:
        """
        Set the orientation of the given atom(s).

        :param rotation: Rotation to assign to the atom(s).
        """
        rot_vec = rotation.rotation_vector
        axis = normalized(rot_vec)
        angle = magnitude(rot_vec) * (radian >> degree)
        self._set(f"quat {axis[0]} {axis[1]} {axis[2]} {angle}")

    def set_shape(self, *, shape: Vector3) -> None:
        """
        Set the shape of the given atom(s).

        :param shape: Shape in LAMMPS units.
        """
        self._set(f"shape {shape[0]} {shape[1]} {shape[2]}")

    def set_peratom_mass(self, *, mass: float) -> None:
        """
        Set the mass of the given atom(s).

        :param mass: Mass in LAMMPS units.
        """
        self._set(f"mass {mass}")

    def set_angular_momentum(self, *, angular_momentum: Vector3) -> None:
        """
        Set the angular momentum of the given atom(s).

        :param angular_momentum: Angular momentum in LAMMPS units.
        """
        self._set(
            f"angmom {angular_momentum[0]} {angular_momentum[1]} {angular_momentum[2]}"
        )


class SingleAtomReference(AtomSetReference):
    """Reference to a single atom by its ID."""

    def __init__(self, simulation: LAMMPSSimulation, atom_id: int):
        super().__init__(simulation)
        self.__id = atom_id

    def _selection_str(self) -> str:
        return f"{self.__id}"

    def _type_str(self) -> str:
        return "atom"


class AtomRangeReference(AtomSetReference):
    """Reference to a range of atoms between two IDs."""

    def __init__(self, simulation: LAMMPSSimulation, atom_id_range: slice):
        super().__init__(simulation)
        self.__range = atom_id_range

    def _selection_str(self) -> str:
        return f"{_slice_to_lammps(self.__range)}"

    def _type_str(self) -> str:
        return "atom"


class _TypeReferenceMixin:
    """Mixin for methods that are valid for atom type references."""

    def set_mass(self: AtomSetReference, mass: float) -> None:  # type: ignore[misc]
        """
        Set the per-type mass.

        :param mass: Mass in LAMMPS units.
        """
        self._simulation.command(f"mass {self._selection_str()} {mass}")


class SingleTypeReference(AtomSetReference, _TypeReferenceMixin):
    """Reference to a single atom type and all atoms with that type."""

    def __init__(self, simulation: LAMMPSSimulation, atom_type: int):
        super().__init__(simulation)
        self.__id = atom_type

    def _selection_str(self) -> str:
        return f"{self.__id}"

    def _type_str(self) -> str:
        return "type"


class TypeRangeReference(AtomSetReference, _TypeReferenceMixin):
    """Reference to a range of atom types and all atoms with those types."""

    def __init__(self, simulation: LAMMPSSimulation, atom_type_range: slice):
        super().__init__(simulation)
        self.__range = atom_type_range

    def _selection_str(self) -> str:
        return f"{_slice_to_lammps(self.__range)}"

    def _type_str(self) -> str:
        return "type"


class SingleMolReference(AtomSetReference):
    """Reference to a single molecule ID and all atoms with that type."""

    def __init__(self, simulation: LAMMPSSimulation, mol_id: int):
        super().__init__(simulation)
        self.__id = mol_id

    def _selection_str(self) -> str:
        return f"{self.__id}"

    def _type_str(self) -> str:
        return "mol"


class MolRangeReference(AtomSetReference):
    """Reference to a range of molecule IDs and all atoms with those types."""

    def __init__(self, simulation: LAMMPSSimulation, mol_id_range: slice):
        super().__init__(simulation)
        self.__range = mol_id_range

    def _selection_str(self) -> str:
        return f"mol {_slice_to_lammps(self.__range)}"

    def _type_str(self) -> str:
        return "mol"


class GroupReference(AtomSetReference):
    """Reference to a group and all atoms within that group."""

    def __init__(self, simulation: LAMMPSSimulation, group_id: int):
        super().__init__(simulation)
        self.__id = group_id

    def _selection_str(self) -> str:
        return f"group {self.__id}"

    def _type_str(self) -> str:
        return "group"
