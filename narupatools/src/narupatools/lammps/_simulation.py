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
from typing import Dict, List, Literal, Optional, TypeVar, Union

import numpy as np
import numpy.typing as npt
from infinite_sets import InfiniteSet
from narupa.trajectory import FrameData

import narupatools.lammps.atom_properties as PROPERTIES
import narupatools.lammps.computes as COMPUTES
import narupatools.lammps.globals as GLOBALS
import narupatools.lammps.settings as SETTINGS
from narupatools.core.units import UnitsNarupa
from narupatools.frame import (
    BondPairs,
    NarupaFrame,
    ParticleCharges,
    ParticleCount,
    ParticleElements,
    ParticleForces,
    ParticleMasses,
    ParticlePositions,
    ParticleResidues,
    ParticleVelocities,
    PotentialEnergy,
)
from narupatools.frame.utils import mass_to_element
from narupatools.lammps._units import get_unit_system
from narupatools.physics.typing import Vector3

from ._constants import VariableDimension, VariableType
from ._wrapper import Extractable, LAMMPSWrapper
from .exceptions import UnknownAtomPropertyError, catch_lammps_warnings_and_exceptions
from .region import Region, RegionSpecification

_TReturnType = TypeVar("_TReturnType")


class LAMMPSSimulation:
    """Wrapper around a LAMMPS simulation."""

    def __init__(self, lammps: LAMMPSWrapper):
        self.__lammps = lammps

        unit_system = get_unit_system(self.extract(GLOBALS.UnitStyle))
        self._lammps_to_narupa = unit_system >> UnitsNarupa
        self._narupa_to_lammps = UnitsNarupa >> unit_system
        self.indexing = LAMMPSIndexing(self)
        self.bond_compute: Optional[COMPUTES.ComputeReference] = None

        self._needs_pre_run = True

        self.mass_compute = self.create_atom_compute(
            id="mass", properties=["mass"], type=VariableType.DOUBLE
        )

        self._kinetic_energy_compute = COMPUTES.KineticEnergy.create(
            self.__lammps, id="ke"
        )

        self._potential_energy_compute = COMPUTES.ComputeGlobalReference.create(
            self.__lammps, id="thermo_pe", dimension=VariableDimension.SCALAR
        )

        self._temperature_compute = COMPUTES.ComputeGlobalReference.create(
            self.__lammps, id="thermo_temp", dimension=VariableDimension.SCALAR
        )
        self._pressure_compute = COMPUTES.ComputeGlobalReference.create(
            self.__lammps, id="thermo_press", dimension=VariableDimension.SCALAR
        )

    @classmethod
    def create_new(
        cls,
        units: Literal["lj", "real", "metal", "si", "cgs", "electron", "micro", "nano"],
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
            lammps.command(f"units {units}")
        return cls(lammps)

    @classmethod
    def from_file(cls, filename: str) -> LAMMPSSimulation:
        """
        Load LAMMPS simulation from a file.

        This automatically runs the 'atom_modify map yes' command required to use atom
        IDs and hence allow certain operations such as setting positions.

        :param filename: Filename of LAMMPS input.
        :return: LAMMPS simulation based on file.
        """
        lammps = LAMMPSWrapper()
        with catch_lammps_warnings_and_exceptions():
            lammps.command("atom_modify map yes")
            lammps.file(filename)
            lammps.command("run 0")
        simulation = cls(lammps)
        return simulation

    def __len__(self) -> int:
        return self.extract(SETTINGS.NumberAllAtoms)

    def gather_atoms(
        self, property: PROPERTIES.AtomProperty[_TReturnType]
    ) -> _TReturnType:
        """
        Gather an atom property from across all processors, and order by Atom ID.

        :param property: The atom property to gather.
        :return: The atom data as a numpy array.
        """
        return self.__lammps.gather_atoms(  # type: ignore[return-value]
            property.key, property.components
        )

    def scatter_atoms(
        self, property: PROPERTIES.AtomProperty, value: np.ndarray
    ) -> None:
        """
        Scatter an atom property across all processors.

        :param property: The atom property to scatter.
        :param value: The value to distribute.
        """
        self.__lammps.scatter_atoms(
            property.key, property.type, property.components, value
        )

    def run(self, steps: int = 1) -> None:
        """
        Run the simulation for a given number of steps.

        :param steps: Number of steps to run.
        """
        pre = "yes" if self._needs_pre_run else "no"
        self.__lammps.command(f"run {steps} pre {pre} post no")
        self._needs_pre_run = False

    @property
    def potential_energy(self) -> float:
        """
        Potential energy of the system in kilojoules per mole.

        The potential energy is computed using the 'thermo_pe' compute which is added to
        all LAMMPS simulations automatically. This uses the 'pe' compute style which
        includes contributions from pairs, bonds, angles, dihedrals, impropers, kspace
        and various fixes.
        """
        return self._potential_energy_compute.extract() * self._lammps_to_narupa.energy

    @property
    def kinetic_energy(self) -> float:
        """Kinetic energy of the system in kilojoules per mole."""
        return self._kinetic_energy_compute.extract() * self._lammps_to_narupa.energy

    @property
    def temperature(self) -> float:
        r"""
        Temperature of the system in kelvin.

        The temperature is computed using the 'thermo_temp' compute which is added to
        all LAMMPS simulations automatically. This uses the 'temp' compute style which
        uses the equation:

        .. math:: KE = \frac{D}{2} N k T

        where :math:`KE` is the kinetic energy, :math:`D` is the dimension of the
        system, :math:`N` is the number of atoms and :math:`k` is the Boltzmann
        constant.

        This removes degrees of freedom due to other fixes that limit molecular motion.
        """
        return self._temperature_compute.extract() * self._lammps_to_narupa.temperature

    @property
    def pressure(self) -> float:
        """Pressure of the system in kilojoules per nanometers cubed."""
        return self._pressure_compute.extract() * self._lammps_to_narupa.pressure

    @property
    def forces(self) -> npt.NDArray[np.float64]:
        """Forces on each atom in kilojoules per mole per nanometer."""
        return self.gather_atoms(PROPERTIES.Force) * self._lammps_to_narupa.force.value

    @property
    def atom_ids(self) -> npt.NDArray[np.int64]:
        """Atom ID of each atom as defined by LAMMPS."""
        return self.gather_atoms(PROPERTIES.AtomID)

    @property
    def charges(self) -> npt.NDArray[np.float64]:
        """Charge of each atom in elementary charges."""
        try:
            return (
                self.gather_atoms(PROPERTIES.Charge)
                * self._lammps_to_narupa.charge.value
            )
        except UnknownAtomPropertyError:
            raise AttributeError

    @property
    def masses(self) -> npt.NDArray[np.float64]:
        """Masses of each atom in daltons."""
        return self.mass_compute.gather()

    @property
    def positions(self) -> npt.NDArray[np.float64]:
        """Positions of each atom in nanometers."""
        return (
            self.gather_atoms(PROPERTIES.Position) * self._lammps_to_narupa.length.value
        )

    @positions.setter
    def positions(self, positions: np.ndarray) -> None:
        self.scatter_atoms(
            PROPERTIES.Position, positions * self._narupa_to_lammps.length
        )

    @property
    def velocities(self) -> npt.NDArray[np.float64]:
        """Velocities of each atom in nanometers per picoseconds."""
        return self.gather_atoms(PROPERTIES.Velocity) * self._lammps_to_narupa.velocity

    @velocities.setter
    def velocities(self, velocities: np.ndarray) -> None:
        self.scatter_atoms(
            PROPERTIES.Velocity, velocities * self._narupa_to_lammps.velocity
        )

    def add_imd_force(self) -> None:
        """Define an IMD force for the LAMMPS simulation."""
        self.command("fix imdprop all property/atom d_imdfx d_imdfy d_imdfz")
        self.command("set atom * d_imdfx 0 d_imdfy 0 d_imdfz 0")
        self.imd_compute = self.create_atom_compute(
            id="imdf",
            properties=["d_imdfx", "d_imdfy", "d_imdfz"],
            type=VariableType.DOUBLE,
        )
        self.command("variable imdfx atom c_imdf[1]")
        self.command("variable imdfy atom c_imdf[2]")
        self.command("variable imdfz atom c_imdf[3]")
        self.command("fix imd_force all addforce v_imdfx v_imdfy v_imdfz")
        self.command("fix_modify imd_force energy yes")

    def get_imd_forces(self) -> npt.NDArray[np.float64]:
        """Get the IMD forces on each atom in kilojoules per mole per angstrom."""
        if self.imd_compute is None:
            raise ValueError("IMD Force not defined.")
        return self.imd_compute.gather() * self._lammps_to_narupa.force

    def set_imd_force(self, index: int, force: Vector3) -> None:
        """
        Set the IMD force on a specific atom.

        :param index: Index of the particle to set.
        :param force: IMD force to apply, in kilojoules per mole per nanometer.
        """
        self.indexing.recompute()
        id = self.indexing.ordered_to_atom_id(index)
        force = force * self._narupa_to_lammps.force
        self.command(
            f"set atom {id} d_imdfx {force[0]} d_imdfy {force[1]} d_imdfz {force[2]}"
        )

    def clear_imd_force(self, index: int) -> None:
        """
        Clear the IMD force on a specific atom.

        :param index: Index of the particle to clear.
        """
        self.indexing.recompute()
        id = self.indexing.ordered_to_atom_id(index)
        self.command(f"set atom {id} d_imdfx 0 d_imdfy 0 d_imdfz 0")

    def command(self, command: str) -> None:
        """
        Run an arbitrary LAMMPS command.

        Running it through this interface means that all previously cached data such as
        positions etc. are cleared. If you are sure running the command will not affect
        any data, set clear_cache to False.

        :param command: LAMMPS command to run.
        :param clear_cache: Should all cached information be cleared.
        """
        self.__lammps.command(command)
        self._needs_pre_run = True

    def get_frame(
        self, *, fields: InfiniteSet[str], existing: Optional[FrameData] = None
    ) -> FrameData:
        """
        Create a Narupa FrameData of the LAMMPS simulation.

        :param fields: Fields to populate of the frame data.
        :param existing: Preexisting frame data.
        :return: FrameData with the appropriate fields.
        """
        if existing is None:
            frame: FrameData = NarupaFrame()
        else:
            frame = existing

        if ParticleCount.key in fields:
            ParticleCount.set(frame, len(self))
        if ParticlePositions.key in fields:
            ParticlePositions.set(frame, self.positions)
        if ParticleVelocities.key in fields:
            ParticleVelocities.set(frame, self.velocities)
        if ParticleMasses.key in fields:
            ParticleMasses.set(frame, self.masses)
        if ParticleElements.key in fields:
            elements = np.vectorize(mass_to_element)(self.masses)
            if np.all(elements is not None):
                ParticleElements.set(frame, elements)
        if ParticleCharges.key in fields:
            with contextlib.suppress(AttributeError):
                ParticleCharges.set(frame, self.charges)
        if ParticleForces.key in fields:
            ParticleForces.set(frame, self.forces)
        if PotentialEnergy.key in fields:
            PotentialEnergy.set(frame, self.potential_energy)

        if BondPairs.key in fields and self.extract(
            SETTINGS.AtomStylesIncludesMolecularTopology
        ):
            self.indexing.recompute()
            if self.bond_compute is None:
                self.bond_compute = self.create_local_compute(
                    id="bonds",
                    properties=["btype", "batom1", "batom2"],
                )
            bonds = self.bond_compute.extract()
            BondPairs.set(frame, self.indexing.atom_id_to_ordered(bonds[:, 1:]))

        if ParticleResidues.key in fields and self.extract(
            SETTINGS.AtomStylesIncludesMolecularTopology
        ):
            self.indexing.recompute()
            mol_ids = self.gather_atoms(PROPERTIES.MoleculeID)

            unique_mol_ids = np.unique(mol_ids)
            mol_id_to_index = {
                mol_id: mol_index
                for mol_index, mol_id in zip(np.argsort(unique_mol_ids), unique_mol_ids)
            }
            ParticleResidues.set(frame, np.vectorize(mol_id_to_index.get)(mol_ids))

        return frame

    def create_local_compute(
        self,
        *,
        id: str,
        properties: List[str],
    ) -> COMPUTES.ComputeReference[npt.NDArray[np.float64]]:
        """Create a local compute."""
        self.__lammps.command(f"compute {id} all property/local {' '.join(properties)}")
        return COMPUTES.ComputeLocalReference(
            self.__lammps,
            id=id,
        )

    def create_atom_compute(
        self,
        *,
        id: str,
        properties: List[str],
        type: Literal[VariableType.INTEGER, VariableType.DOUBLE],
    ) -> COMPUTES.ComputeAtomReference[npt.NDArray[np.float64]]:
        """
        Create a per-atom compute.

        :param id: Compute ID to do.
        :param properties: Atom properties to compute.
        :param type: Variable type to return.
        :return: Reference to the compute.
        """
        self.command(f"compute {id} all property/atom {' '.join(properties)}")
        return COMPUTES.ComputeAtomReference(
            self.__lammps,
            id=id,
            type=type,
            count=len(properties),
        )

    @property
    def timestep(self) -> float:
        """Timestep of the simulation in picoseconds."""
        return self.extract(GLOBALS.TimestepLength) * self._lammps_to_narupa.time

    @timestep.setter
    def timestep(self, value: float) -> None:
        self.command(f"timestep {value * self._narupa_to_lammps.time}")
        self._timestep = None

    def create_atom(self, *, type: int, position: Vector3) -> None:
        """
        Insert an atom into the simulation.

        :param type: Type of the atom to insert.
        :param position: Position of the atom in nanometers.
        """
        position = position * self._narupa_to_lammps.length
        self.command(
            f"create_atoms {type} single {position[0]} {position[1]} {position[2]}"
        )

    def create_region(
        self, region: RegionSpecification, *, id: Optional[str] = None
    ) -> Region:
        """
        Create a region in the simulation.

        :param region: Specification to the region.
        :param id: Optional id for the region.
        :return: Region object.
        """
        if id is None:
            id = self.__lammps.regions.generate_id()
        self.command(
            f"region {id} {region.style} {region.args(self._narupa_to_lammps)}"
        )
        return Region(self, id, region)

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
        self.command(f"create_box {n_types} {region.id}")

    def set_mass(self, *, type: Union[int, slice], mass: float) -> None:
        """
        Set the mass of a specified atom type.

        :param type: Atom type.
        :param mass: Mass in daltons.
        """
        self.command(f"mass {type} {mass * self._narupa_to_lammps.mass}")

    def setup_langevin(self, *, temperature: float, friction: float, seed: int) -> None:
        """
        Setup a Langevin thermostat.

        :param temperature: Temperature in kelvin
        :param friction: Friction of the thermostat in inverse picoseconds.
        :param seed: Seed to use for random number generation.
        """
        temperature *= self._narupa_to_lammps.temperature
        friction /= self._narupa_to_lammps.time
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
        self.__simulation = simulation
        self.__atom_ids = simulation.create_atom_compute(
            id="atomids", properties=["id"], type=VariableType.INTEGER
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
