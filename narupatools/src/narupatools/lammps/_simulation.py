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
import re
import warnings
from contextlib import contextmanager
from ctypes import Array, c_double
from typing import Any, Dict, Generator, List, Literal, Optional, Union, overload

import numpy as np
from infinite_sets import InfiniteSet
from lammps import PyLammps
from MDAnalysis import Universe
from narupa.trajectory import FrameData

from narupatools.core.units import UnitsNarupa
from narupatools.frame.fields import (
    ParticleCharges,
    ParticleCount,
    ParticleForces,
    ParticleMasses,
    ParticlePositions,
    ParticleVelocities,
    PotentialEnergy,
)
from narupatools.lammps._units import get_unit_system
from narupatools.mdanalysis import mdanalysis_universe_to_frame
from narupatools.physics.typing import Vector3

from ._constants import PropertyType, VariableStyle, VariableType
from ._region import Region, RegionSpecification
from ._warnings import LAMMPSWarning
from .exceptions import (
    AtomIDsNotDefinedError,
    CannotOpenFileError,
    ComputeNotFoundError,
    IllegalCommandError,
    InvalidComputeSpecificationError,
    LAMMPSError,
    MissingInputScriptError,
    UnknownAtomPropertyError,
    UnknownCommandError,
    UnknownPropertyNameError,
    UnrecognizedStyleError,
)
from .output_capture import OutputCapture


class LAMMPSSimulation:
    """Wrapper around a LAMMPS simulation."""

    def __init__(self, lammps: PyLammps, universe: Universe = None):
        self.__lammps = lammps
        self.__lammps.enable_cmd_history = True
        self._universe = universe
        self._index_to_id: Dict[int, int] = {}
        unit_system = get_unit_system(self.__lammps.system.units)
        self._lammps_to_narupa = unit_system >> UnitsNarupa
        self._narupa_to_lammps = UnitsNarupa >> unit_system

        self._potential_energy: Optional[float] = None
        self._kinetic_energy: Optional[float] = None
        self._forces: Optional[np.ndarray] = None
        self._charges: Optional[np.ndarray] = None
        self._positions: Optional[np.ndarray] = None
        self._velocities: Optional[np.ndarray] = None
        self._atom_ids: Optional[np.ndarray] = None
        self._masses: Optional[np.ndarray] = None
        self._size: Optional[int] = None
        self._temperature: Optional[float] = None
        self._pressure: Optional[float] = None
        self._needs_pre_run = True
        self._timestep: Optional[float] = None

        self.command("compute mass all property/atom mass")
        self.command("compute thermo_ke all ke")

        self._clear_cache()

    @property
    def command_history(self) -> List[str]:
        """List of commands run on the simulation."""
        return self.__lammps._cmd_history

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
        lammps = PyLammps()
        with catch_lammps_warnings_and_exceptions():
            lammps.atom_modify("map yes")
            lammps.units(units)
        return cls(lammps)

    @classmethod
    def from_file(cls, filename: str, data_filename: str) -> LAMMPSSimulation:
        """
        Load LAMMPS simulation from a file.

        This automatically runs the 'atom_modify map yes' command required to use atom
        IDs and hence allow certain operations such as setting positions.

        :param filename: Filename of LAMMPS input.
        :param data_filename: Filename of LAMMPS data file, used to generate MDAnalysis
                              universe.
        :return: LAMMPS simulation based on file.
        """
        lammps = PyLammps()
        with catch_lammps_warnings_and_exceptions():
            lammps.atom_modify("map yes")
            lammps.file(filename)
            lammps.command("run 0")
        universe = Universe(data_filename, format="DATA")
        simulation = cls(lammps, universe)
        return simulation

    def eval(self, expression: str) -> Any:
        """
        Evaluate an arbitrary expression.

        :param expression: Expression to be evaluated using the eval command.
        """
        return self.__lammps.eval(expression)

    @property
    def universe(self) -> Universe:
        """An MDAnalysis Universe containing topology information."""
        if self._universe is None:
            raise AttributeError("Simulation does not have MDAnalysis universe")
        return self._universe

    def _clear_cache(self) -> None:
        self._positions = None
        self._velocities = None
        self._forces = None
        self._charges = None
        self._masses = None
        self._atom_ids = None
        self._potential_energy = None
        self._needs_pre_run = True
        self._temperature = None
        self._size = None

    def __len__(self) -> int:
        if self._size is None:
            self._size = self.__lammps.lmp.get_natoms()
        return self._size

    @overload
    def extract_compute(
        self, key: str, style: VariableStyle, type: Literal[VariableType.SCALAR]
    ) -> float:
        ...

    @overload
    def extract_compute(
        self, key: str, style: VariableStyle, type: Literal[VariableType.ARRAY]
    ) -> np.ndarray:
        ...

    @overload
    def extract_compute(
        self, key: str, style: VariableStyle, type: Literal[VariableType.VECTOR]
    ) -> np.ndarray:
        ...

    def extract_compute(
        self, key: str, style: VariableStyle, type: VariableType
    ) -> Union[float, np.ndarray]:
        """
        Extract value of a compute.

        :param key: ID of the compute.
        :param style: Style of the variable.
        :param type: Type of the variable.
        :raises ComputeNotFoundError: Key does not match any known computes.
        :raises InvalidComputeSpecificationError: Compute exists but doesn't match style
                                                  and type provided.
        :return: Value of the compute.
        """
        value = self.__lammps.lmp.numpy.extract_compute(key, style, type)
        if value is None:
            computes = self.__lammps.computes
            if not any(compute["name"] == key for compute in computes):
                raise ComputeNotFoundError(key)
            else:
                raise InvalidComputeSpecificationError(key, style, type)
        return value

    def gather_atoms(self, key: str, type: PropertyType, dimension: int) -> np.ndarray:
        """
        Gather a per-atom property from all the processors.

        :param key: Id of the property.
        :param type: Type of the variable.
        :param dimension: Dimension of the variable.
        :raises UnknownAtomPropertyError: Property was not found in the simulation.
        :return: Value of the atom property.
        """
        dtype = float if type == PropertyType.DOUBLE else int
        try:
            with catch_lammps_warnings_and_exceptions():
                natoms = self.__lammps.lmp.get_natoms()
                data = ((dimension * natoms) * type.get_ctype())()

                self.__lammps.lmp.lib.lammps_gather_atoms(  # type: ignore[attr-defined]
                    self.__lammps.lmp.lmp, key.encode(), type, dimension, data  # type: ignore[attr-defined]
                )

            if dimension == 1:
                return np.array(data, dtype=dtype).reshape((-1,))
            else:
                return np.array(data, dtype=dtype).reshape((-1, dimension))
        except UnknownPropertyNameError:
            # Reraise as a different error so we know what key caused the error.
            raise UnknownAtomPropertyError(key)

    def _scatter_atoms(
        self, name: str, type: PropertyType, dimensions: int, value: np.ndarray
    ) -> None:
        n_atoms = self.__lammps.lmp.get_natoms()
        with catch_lammps_warnings_and_exceptions():
            self.__lammps.lmp.scatter_atoms(
                name, type, dimensions, _to_ctypes(value, n_atoms)
            )
        self._clear_cache()

    def run(self, steps: int = 1) -> None:
        """
        Run the simulation for a given number of steps.

        :param steps: Number of steps to run.
        """
        pre = "yes" if self._needs_pre_run else "no"
        self.command(f"run {steps} pre {pre} post no")
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
        if self._potential_energy is None:
            self._potential_energy = (
                self.extract_compute(
                    "thermo_pe", VariableStyle.GLOBAL, VariableType.SCALAR
                )
                * self._lammps_to_narupa.energy
            )
        return self._potential_energy

    @property
    def kinetic_energy(self) -> float:
        """Kinetic energy of the system in kilojoules per mole."""
        if self._kinetic_energy is None:
            self._kinetic_energy = (
                self.extract_compute(
                    "thermo_ke", VariableStyle.GLOBAL, VariableType.SCALAR
                )
                * self._lammps_to_narupa.energy
            )
        return self._kinetic_energy

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
        if self._temperature is None:
            self._temperature = (
                self.extract_compute(
                    "thermo_temp", VariableStyle.GLOBAL, VariableType.SCALAR
                )
                * self._lammps_to_narupa.temperature
            )
        return self._temperature

    @property
    def pressure(self) -> float:
        """Pressure of the system in kilojoules per nanometers cubed."""
        if self._pressure is None:
            self._pressure = (
                self.extract_compute(
                    "thermo_press", VariableStyle.GLOBAL, VariableType.SCALAR
                )
                * self._lammps_to_narupa.pressure
            )
        return self._pressure

    @property
    def forces(self) -> np.ndarray:
        """Forces on each atom in kilojoules per mole per nanometer."""
        if self._forces is None:
            self._forces = (
                self.gather_atoms("f", PropertyType.DOUBLE, 3)
                * self._lammps_to_narupa.force
            )
        return self._forces

    @property
    def atom_ids(self) -> np.ndarray:
        """ID of each atom."""
        if self._atom_ids is None:
            self._atom_ids = self.gather_atoms("id", PropertyType.INT, 1)
        return self._atom_ids

    @property
    def charges(self) -> np.ndarray:
        """Charge of each atom in elementary charges."""
        if self._charges is None:
            try:
                self._charges = (
                    self.gather_atoms("q", PropertyType.DOUBLE, 1)
                    * self._lammps_to_narupa.charge
                )
            except UnknownAtomPropertyError:
                raise AttributeError
        return self._charges

    @property
    def masses(self) -> np.ndarray:
        """Masses of each atom in daltons."""
        if self._masses is None:
            self._masses = (
                self.extract_compute("mass", VariableStyle.ATOM, VariableType.VECTOR)
                * self._lammps_to_narupa.mass
            )
        return self._masses

    @property
    def positions(self) -> np.ndarray:
        """Positions of each atom in nanometers."""
        if self._positions is None:
            self._positions = (
                self.gather_atoms("x", PropertyType.DOUBLE, 3)
                * self._lammps_to_narupa.length
            )
        return self._positions

    @positions.setter
    def positions(self, positions: np.ndarray) -> None:
        self._scatter_atoms(
            "x", PropertyType.DOUBLE, 3, positions * self._narupa_to_lammps.length
        )

    @property
    def velocities(self) -> np.ndarray:
        """Velocities of each atom in nanometers per picoseconds."""
        if self._velocities is None:
            self._velocities = (
                self.gather_atoms("v", PropertyType.DOUBLE, 3)
                * self._lammps_to_narupa.velocity
            )
        return self._velocities

    @velocities.setter
    def velocities(self, velocities: np.ndarray) -> None:
        self._scatter_atoms(
            "v", PropertyType.DOUBLE, 3, velocities * self._narupa_to_lammps.velocity
        )

    def get_imd_forces(self) -> np.ndarray:
        """Get the IMD forces on each atom in kilojoules per mole per angstrom."""
        imd_extracted = self.__lammps.lmp.numpy.extract_compute(
            "imdf", VariableStyle.ATOM, VariableType.ARRAY
        )
        return imd_extracted * self._lammps_to_narupa.force  # type: ignore

    def _generate_index_to_id(self) -> None:
        ids = self.atom_ids
        for index, id in enumerate(ids):
            self._index_to_id[index] = id

    def set_imd_force(self, index: int, force: Vector3) -> None:
        """
        Set the IMD force on a specific atom.

        :param index: Index of the particle to set.
        :param force: IMD force to apply, in kilojoules per mole per nanometer.
        """
        if len(self._index_to_id) == 0:
            self._generate_index_to_id()
        id = self._index_to_id[index]
        force = force * self._narupa_to_lammps.force
        self.command(
            f"set atom {id} d_imdfx {force[0]} d_imdfy {force[1]} d_imdfz {force[2]}"
        )
        self._clear_cache()

    def clear_imd_force(self, index: int) -> None:
        """
        Clear the IMD force on a specific atom.

        :param index: Index of the particle to clear.
        """
        if len(self._index_to_id) == 0:
            self._generate_index_to_id()
        id = self._index_to_id[index]
        self.command(f"set atom {id} d_imdfx 0 d_imdfy 0 d_imdfz 0")
        self._clear_cache()

    def command(self, command: str, clear_cache: bool = True) -> None:
        """
        Run an arbitrary LAMMPS command.

        Running it through this interface means that all previously cached data such as
        positions etc. are cleared. If you are sure running the command will not affect
        any data, set clear_cache to False.

        :param command: LAMMPS command to run.
        :param clear_cache: Should all cached information be cleared.
        """
        with catch_lammps_warnings_and_exceptions():
            self.__lammps.command(command)
        if clear_cache:
            self._clear_cache()

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
            frame = FrameData()
        else:
            frame = existing

        if self._universe is not None:
            frame = mdanalysis_universe_to_frame(
                self._universe, fields=fields, frame=frame
            )
        if ParticleCount.key in fields:
            ParticleCount.set(frame, len(self))
        if ParticlePositions.key in fields:
            ParticlePositions.set(frame, self.positions)
        if ParticleVelocities.key in fields:
            ParticleVelocities.set(frame, self.velocities)
        if ParticleMasses.key in fields:
            ParticleMasses.set(frame, self.masses)
        if ParticleCharges.key in fields:
            with contextlib.suppress(AttributeError):
                ParticleCharges.set(frame, self.charges)
        if ParticleForces.key in fields:
            ParticleForces.set(frame, self.forces)
        if PotentialEnergy.key in fields:
            PotentialEnergy.set(frame, self.potential_energy)
        return frame

    @property
    def lammps_packages(self) -> List[str]:
        """List of installed LAMMPS packages."""
        return self.__lammps.lmp.installed_packages

    @property
    def computes(self) -> Dict[str, Compute]:
        """Dictionary of computes indexed by their names."""
        return {
            compute["name"]: Compute(
                simulation=self,
                group=compute["group"],
                name=compute["name"],
                style=compute["style"],
            )
            for compute in self.__lammps.computes
        }

    @property
    def timestep(self) -> float:
        """Timestep of the simulation in picoseconds."""
        if self._timestep is None:
            self._timestep = self.eval("dt") * self._lammps_to_narupa.time
        return self._timestep

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
            x = 1
            ids = self.region_ids
            while str(x) in ids:
                x += 1
            id = str(x)
        self.command(
            f"region {id} {region.style} {region.args(self._narupa_to_lammps)}"
        )
        return Region(self, id, region)

    @property
    def region_ids(self) -> List[str]:
        """List of region ids currently in use for the simulation."""
        ids = []
        info = self.__lammps.info("regions")
        pattern = re.compile(r"Region\[\s*\d+\]:\s*([\w\-_]+)")
        for line in info:
            if (match := pattern.match(line)) is not None:
                ids.append(match.group(1))
        return ids

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

    def add_imd_force(self) -> None:
        """Define an IMD force for the LAMMPS simulation."""
        self.command("fix imdprop all property/atom d_imdfx d_imdfy d_imdfz")
        self.command("set atom * d_imdfx 0 d_imdfy 0 d_imdfz 0")
        self.command("compute imdf all property/atom d_imdfx d_imdfy d_imdfz")
        self.command("variable imdfx atom c_imdf[1]")
        self.command("variable imdfy atom c_imdf[2]")
        self.command("variable imdfz atom c_imdf[3]")
        self.command("fix imd_force all addforce v_imdfx v_imdfy v_imdfz")
        self.command("fix_modify imd_force energy yes")

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


class Compute:
    """Represents a compute defined in a LAMMPS simulation."""

    def __init__(
        self, *, simulation: LAMMPSSimulation, group: str, name: str, style: str
    ):
        self._simulation = simulation
        self._group = group
        self._name = name
        self._style = style

    def __eq__(self, other: Any) -> bool:
        return (
            isinstance(other, Compute)
            and self._simulation is other._simulation
            and self._group == other._group
            and self._name == other._name
            and self._style == other._style
        )

    def __hash__(self) -> int:
        return hash((self._simulation, self._group, self._name, self._style))


def _to_ctypes(array: np.ndarray, natoms: int) -> Array:
    n3 = 3 * natoms
    x = (n3 * c_double)()
    for i, f in enumerate(array.flat):
        x[i] = f
    return x


ILLEGAL_COMMAND_REGEX = re.compile(r"illegal \w+ command", re.IGNORECASE)
UNKNOWN_COMMAND_REGEX = re.compile(r"unknown command", re.IGNORECASE)
MISSING_INPUT_SCRIPT_REGEX = re.compile(
    r"cannot open input script [\w.]+: No such file or directory", re.IGNORECASE
)
CANNOT_OPEN_FILE_REGEX = re.compile(
    r"cannot open file [\w.]+: No such file or directory", re.IGNORECASE
)
UNRECOGNIZED_STYLE_REGEX = re.compile(r"Unrecognized \w+ style", re.IGNORECASE)


def _handle_error(message: str) -> None:
    if MISSING_INPUT_SCRIPT_REGEX.match(message) is not None:
        raise MissingInputScriptError(message)
    if CANNOT_OPEN_FILE_REGEX.match(message) is not None:
        raise CannotOpenFileError(message)
    if ILLEGAL_COMMAND_REGEX.match(message) is not None:
        raise IllegalCommandError(message)
    if UNKNOWN_COMMAND_REGEX.match(message) is not None:
        raise UnknownCommandError(message)
    if UNRECOGNIZED_STYLE_REGEX.match(message) is not None:
        raise UnrecognizedStyleError(message)
    raise LAMMPSError(message)


@contextmanager
def catch_lammps_warnings_and_exceptions() -> Generator[None, None, None]:
    """
    Capture output and raises logged warnings and errors in a pythonic way.

    Any line starting with 'WARNING: ' will be raised as a LAMMPSWarning, except certain
    warnings which are instead treated as errors.

    Any line starting with 'ERROR:' will be raised as a LAMMPSError, unless a more
    specific error message exists.

    :raises AtomIDsNotDefinedError: Library error raises in gather/scatter atoms.
    :raises LAMMPSError: An ERROR: is logged to the console by LAMMPS.
    """
    with OutputCapture() as o:
        try:
            yield
        except Exception as e:
            if e.args[0].startswith("ERROR on proc 0: "):
                _handle_error(e.args[0][17:])
            if e.args[0].startswith("ERROR: "):
                _handle_error(e.args[0][7:])
            raise LAMMPSError(e.args[0])
        output = o.output
    for line in output.splitlines():
        if line.startswith("WARNING: "):
            warning = line[9:]
            if warning.startswith("Library error in lammps_gather_atoms"):
                raise AtomIDsNotDefinedError(func_name="gather_atoms")
            elif warning.startswith("Library error in lammps_scatter_atoms"):
                raise AtomIDsNotDefinedError(func_name="scatter_atoms")
            else:
                warnings.warn(LAMMPSWarning(warning))
        if line.startswith("ERROR: "):
            error = line[7:]
            _handle_error(error)
