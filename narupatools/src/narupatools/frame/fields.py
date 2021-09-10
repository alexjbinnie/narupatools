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

"""Standard fields and their getters/setters for Narupa frames."""

from abc import ABCMeta, abstractmethod
from typing import Dict, Generic, Iterable, TypeVar, Union

import numpy as np
import numpy.typing as npt
from narupa.trajectory.frame_data import (
    BOND_ORDERS,
    BOND_PAIRS,
    BOX_VECTORS,
    CHAIN_COUNT,
    CHAIN_NAMES,
    KINETIC_ENERGY,
    PARTICLE_COUNT,
    PARTICLE_ELEMENTS,
    PARTICLE_NAMES,
    PARTICLE_POSITIONS,
    PARTICLE_RESIDUES,
    PARTICLE_TYPES,
    POTENTIAL_ENERGY,
    RESIDUE_CHAINS,
    RESIDUE_COUNT,
    RESIDUE_IDS,
    RESIDUE_NAMES,
    FrameData,
)
from typing_extensions import Final

from narupatools.frame._utils import atomic_numbers_to_masses

PARTICLE_MASSES = "particle.masses"
PARTICLE_VELOCITIES = "particle.velocities"
PARTICLE_FORCES = "particle.forces"
BOND_COUNT = "bond.count"
BOND_TYPES = "bond.types"

_TFrom = TypeVar("_TFrom")
_TTo = TypeVar("_TTo")

AssignableToFloatArray = npt.ArrayLike
AssignableToIndexArray = npt.ArrayLike
AssignableToStringArray = Iterable[str]


def _to_n_by_3(value: npt.ArrayLike) -> np.ndarray:
    """
    Convert a generic float array to a 3 by N NumPy array.

    This is a faster implementation that what is currently in Narupa.

    :param value: Object that can be converted to a NumPy array.
    """
    return np.array(value).reshape((-1, 3))


def _to_n_by_2(value: npt.ArrayLike) -> np.ndarray:
    """
    Convert a generic float array to a 2 by N NumPy array.

    This is a faster implementation that what is currently in Narupa.

    :param value: Object that can be converted to a NumPy array.
    """
    return np.array(value).reshape((-1, 2))


def _flatten_array(value: npt.ArrayLike) -> np.ndarray:
    return np.asarray(value).ravel()


_DEFINED_KEYS: Dict[str, "FrameKey"] = {}


def get_frame_key(key: str) -> "FrameKey":
    """
    Get the defined FrameKey for a given string key.

    :param key: Key as it appears in the frame data.
    """
    return _DEFINED_KEYS[key]


_TDefault = TypeVar("_TDefault")


class FrameKey(Generic[_TFrom, _TTo], metaclass=ABCMeta):
    """
    Key that can be placed in a Narupa FrameData, such as positions or potential energy.

    This allows custom keys to be more easily defined, instead of the hardcoded
    shortcuts of the original Narupa FrameData.
    """

    key: Final[str]

    def __init__(self, key: str):
        global _DEFINED_KEYS
        self.key = key
        _DEFINED_KEYS[key] = self

    @abstractmethod
    def set(self, frame_data: FrameData, value: _TFrom) -> None:
        """
        Insert the value into the given frame data for this key.

        :param frame_data: FrameData to modify.
        :param value: Value to insert.
        """
        ...

    @abstractmethod
    def _get(self, frame_data: FrameData, /) -> _TTo:
        ...

    def _calculate(self, frame: FrameData, /) -> _TTo:
        """
        Attempt to calculate the value based on other keys.

        :param frame: Frame to look for other keys:
        :raises KeyError: Cannot calculate value based on other available keys.
        """
        raise KeyError

    def get(self, frame_data: FrameData, /, *, calculate: bool = False) -> _TTo:
        """
        Get the value of the given frame data for this key.

        :param frame_data: FrameData to lookup.
        :param calculate: Calculate the value based on other keys available.
        :raises KeyError: Given key is not available.
        """
        try:
            return self._get(frame_data)
        except KeyError:
            if calculate:
                return self._calculate(frame_data)
        raise KeyError(f"Frame does not contain key {self.key}")

    def get_with_default(
        self, frame_data: FrameData, /, default: _TDefault, *, calculate: bool = False
    ) -> Union[_TTo, _TDefault]:
        """
        Get the value for this key, returning a default value if the key is absent.

        :param frame_data: FrameData to look up key in.
        :param default: Default value to return if absent.
        :param calculate: Calculate the value based on other keys available.
        """
        try:
            return self.get(frame_data, calculate=calculate)
        except KeyError:
            return default


class _FloatArrayKey(FrameKey[AssignableToFloatArray, np.ndarray]):
    def set(self, frame_data: FrameData, value: AssignableToFloatArray) -> None:
        frame_data.set_float_array(self.key, value)

    def _get(self, frame_data: FrameData) -> np.ndarray:
        if self.key not in frame_data.raw.arrays.keys():
            raise KeyError(f"{self.key} not present in FrameData")
        return np.array(frame_data.raw.arrays[self.key].float_values.values)


class _ThreeByNFloatArrayKey(FrameKey[AssignableToFloatArray, np.ndarray]):
    def set(self, frame_data: FrameData, value: AssignableToFloatArray) -> None:
        array = _flatten_array(value)
        if len(array) % 3 > 0:
            raise TypeError(f"Cannot set {self.key} to array not divisible by 3.")
        frame_data.set_float_array(self.key, array)

    def _get(self, frame_data: FrameData) -> np.ndarray:
        if self.key not in frame_data.raw.arrays.keys():
            raise KeyError(f"{self.key} not present in FrameData")
        return _to_n_by_3(frame_data.raw.arrays[self.key].float_values.values)


class _IntegerArrayKey(FrameKey[AssignableToIndexArray, np.ndarray]):
    def set(self, frame_data: FrameData, value: AssignableToIndexArray) -> None:
        frame_data.set_index_array(self.key, value)

    def _get(self, frame_data: FrameData) -> np.ndarray:
        if self.key not in frame_data.raw.arrays.keys():
            raise KeyError(f"{self.key} not present in FrameData")
        return np.array(frame_data.raw.arrays[self.key].index_values.values)


class _TwoByNIntegerArrayKey(FrameKey[AssignableToIndexArray, np.ndarray]):
    def set(self, frame_data: FrameData, value: AssignableToIndexArray) -> None:
        array = _flatten_array(value)
        if len(array) % 2 > 0:
            raise TypeError(f"Cannot set {self.key} to array not divisible by 2.")
        frame_data.set_index_array(self.key, array)

    def _get(self, frame_data: FrameData) -> np.ndarray:
        if self.key not in frame_data.raw.arrays.keys():
            raise KeyError(f"{self.key} not present in FrameData")
        return _to_n_by_2(frame_data.raw.arrays[self.key].index_values.values)


class _StringArrayKey(FrameKey[AssignableToStringArray, np.ndarray]):
    def set(self, frame_data: FrameData, value: AssignableToStringArray) -> None:
        frame_data.set_string_array(self.key, value)

    def _get(self, frame_data: FrameData) -> np.ndarray:
        if self.key not in frame_data.raw.arrays.keys():
            raise KeyError(f"{self.key} not present in FrameData")
        return np.array(frame_data.raw.arrays[self.key].string_values.values)


class _IntegerKey(FrameKey[int, int]):
    def set(self, frame_data: FrameData, value: int) -> None:
        frame_data.raw.values[self.key].number_value = value

    def _get(self, frame_data: FrameData) -> int:
        if self.key not in frame_data.raw.values.keys():
            raise KeyError(f"{self.key} not present in FrameData")
        return int(frame_data.raw.values[self.key].number_value)


class _FloatKey(FrameKey[float, float]):
    def set(self, frame_data: FrameData, value: float) -> None:
        frame_data.raw.values[self.key].number_value = value

    def _get(self, frame_data: FrameData) -> float:
        if self.key not in frame_data.raw.values.keys():
            raise KeyError(f"{self.key} not present in FrameData")
        return float(frame_data.raw.values[self.key].number_value)


class _ParticleMassesKey(_FloatArrayKey):
    def _calculate(self, frame_data: FrameData, /) -> np.ndarray:
        elements = ParticleElements.get(frame_data)
        return atomic_numbers_to_masses(elements)


ParticlePositions = _ThreeByNFloatArrayKey(PARTICLE_POSITIONS)
"""
Array of particle positions in nanometers, as a NumPy N by 3 array of floats.
"""

ParticleVelocities = _ThreeByNFloatArrayKey(PARTICLE_VELOCITIES)
"""
Array of particle velocities in nanometers per picosecond, as a NumPy N by 3 array of
floats.
"""

ParticleForces = _ThreeByNFloatArrayKey(PARTICLE_FORCES)
"""
Array of particle forces in kilojoule per mole per nanometer, as a NumPy N by 3 array of
floats.
"""

ParticleElements = _IntegerArrayKey(PARTICLE_ELEMENTS)
"""
Array of particle elements, as a NumPy array of atomic number.
"""

ParticleResidues = _IntegerArrayKey(PARTICLE_RESIDUES)
"""
Array of particle residue indices, as a NumPy array of zero-based indices.
"""

ParticleNames = _StringArrayKey(PARTICLE_NAMES)
"""
Array of particle names, as a NumPy array of strings.
"""

ParticleTypes = _StringArrayKey(PARTICLE_TYPES)
"""
Array of particle types, as a NumPy array of strings.
"""

ParticleMasses = _ParticleMassesKey(PARTICLE_MASSES)
"""
Array of particle masses in atomic mass units, as a NumPy array of floats.
"""

ParticleCharges = _FloatArrayKey("particle.charges")
"""
Array of particle charges in proton charges, as a NumPy array of floats.
"""

ParticleCount = _IntegerKey(PARTICLE_COUNT)
"""
Number of particles in the system, as an integer.
"""

PotentialEnergy = _FloatKey(POTENTIAL_ENERGY)
"""
Potential energy of the system in kilojoules per mole, as a float.
"""

KineticEnergy = _FloatKey(KINETIC_ENERGY)
"""
Kinetic energy of the system in kilojoules per mole, as a float.
"""

SimulationElapsedTime = _FloatKey("simulation.elapsed_time")
"""
Elapsed time in the simulation in picoseconds, as a float.
"""

SimulationTotalTime = _FloatKey("simulation.total_time")
"""
Total time in the simulation in picoseconds, as a float.
"""

SimulationElapsedSteps = _IntegerKey("simulation.elapsed_steps")
"""
Elapsed steps in the simulation, as a integer.
"""

SimulationTotalSteps = _IntegerKey("simulation.total_steps")
"""
Total steps in the simulation, as a integer.
"""

BondPairs = _TwoByNIntegerArrayKey(BOND_PAIRS)
"""
Array of bonds as pairs of particle indices, as a N by 2 NumPy array of integers.
"""

BondOrders = _IntegerArrayKey(BOND_ORDERS)
"""
Array of bond orders, as a NumPy array of integers.
"""

BondTypes = _StringArrayKey(BOND_TYPES)
"""
Array of bond types, as a NumPy array of strings.
"""

BondCount = _IntegerKey(BOND_COUNT)
"""
Number of bonds in the system, as an integer.
"""

ResidueNames = _StringArrayKey(RESIDUE_NAMES)
"""
Array of residue names, as a NumPy array of strings.
"""

ResidueIds = _StringArrayKey(RESIDUE_IDS)
"""
Array of residue ids, as a NumPy array of strings.
"""

ResidueChains = _IntegerArrayKey(RESIDUE_CHAINS)
"""
Array of residue chain indices, as a NumPy array of integers.
"""

ResidueCount = _IntegerKey(RESIDUE_COUNT)
"""
Number of residues in the system, as an integer.
"""

ChainNames = _StringArrayKey(CHAIN_NAMES)
"""
Array of chain names, as a NumPy array of strings.
"""

ChainCount = _IntegerKey(CHAIN_COUNT)
"""
Number of chains in the system, as an integer.
"""

BoxVectors = _ThreeByNFloatArrayKey(BOX_VECTORS)
"""
Array of box vectors in nanometers indicating the simulation box, as a NumPy 3 by 3 or
4 x 3 array of floats.

If a 3 x 3 matrix is provided, the columns give the directions and magnitudes of the
three box vectors. If a fourth column is present, this indicates the translation of the
box.
"""

DYNAMIC_FIELDS = frozenset(
    {
        ParticlePositions.key,
        PotentialEnergy.key,
        KineticEnergy.key,
        ParticleVelocities.key,
        ParticleForces.key,
        SimulationTotalSteps.key,
        SimulationTotalTime.key,
        SimulationElapsedTime.key,
        SimulationElapsedSteps.key,
    }
)
"""Set of fields that are by default considered to change when a simulation is run."""
