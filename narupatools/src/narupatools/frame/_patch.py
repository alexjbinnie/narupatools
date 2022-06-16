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

"""Patch to FrameData that allows copy() to work with empty arrays."""

from typing import Any, Dict, Generator, ItemsView, KeysView, Union

import numpy as np
from infinite_sets import InfiniteSet, everything
from narupa.trajectory import FrameData
from narupa.utilities.protobuf_utilities import value_to_object

from narupatools.physics.typing import ScalarArray, Vector3Array
from narupatools.state.typing import Serializable
from narupatools.util import monkeypatch

from ._properties import DynamicStructureMethods
from .fields import (
    BondCount,
    BondOrders,
    BondPairs,
    BondTypes,
    BoxVectors,
    ChainCount,
    ChainNames,
    FrameKey,
    KineticEnergy,
    ParticleAngularMomenta,
    ParticleCount,
    ParticleElements,
    ParticleForces,
    ParticleMasses,
    ParticleNames,
    ParticlePositions,
    ParticlePrincipalMoments,
    ParticleResidues,
    ParticleRotations,
    ParticleTypes,
    ParticleVelocities,
    PotentialEnergy,
    ResidueChains,
    ResidueCount,
    ResidueIds,
    ResidueNames,
    get_frame_key,
)


@monkeypatch(FrameData)
class _PatchedFrameData(DynamicStructureMethods, FrameData):
    bond_pairs = BondPairs  # type: ignore
    bond_orders = BondOrders  # type: ignore
    bond_types = BondTypes
    bond_count = BondCount

    particle_positions = ParticlePositions  # type: ignore
    particle_rotations = ParticleRotations
    particle_elements = ParticleElements  # type: ignore
    particle_types = ParticleTypes  # type: ignore
    particle_names = ParticleNames  # type: ignore
    particle_residues = ParticleResidues  # type: ignore
    particle_masses = ParticleMasses
    particle_velocities = ParticleVelocities
    particle_forces = ParticleForces
    particle_count = ParticleCount  # type: ignore
    particle_angular_momenta = ParticleAngularMomenta
    particle_principal_moments = ParticlePrincipalMoments

    residue_names = ResidueNames  # type: ignore
    residue_uds = ResidueIds
    residue_chains = ResidueChains  # type: ignore
    residue_count = ResidueCount  # type: ignore

    chain_names = ChainNames  # type: ignore
    chain_count = ChainCount  # type: ignore

    kinetic_energy = KineticEnergy  # type: ignore
    potential_energy = PotentialEnergy  # type: ignore
    box_vectors = BoxVectors  # type: ignore

    def copy(self, fields: InfiniteSet[str] = everything()) -> FrameData:
        frame = FrameData()
        for key, value in self.raw.arrays.items():
            if key in fields:
                frame.raw.arrays[key].CopyFrom(value)
        for key, value in self.raw.values.items():
            if key in fields:
                frame.raw.values[key].CopyFrom(value)
        return frame

    def __repr__(self) -> str:
        rep = "<FrameData"

        for key, value in self.items():
            rep += f" {key}={_print_value(value)}"

        rep += ">"

        return rep

    def keys(self) -> KeysView[str]:
        """Iterate over the keys in the Frame."""
        return KeysView(self)  # type: ignore

    def items(self) -> ItemsView[str, Any]:
        """Iterate over the keys and values of the Frame."""
        return ItemsView(self)  # type: ignore

    def __iter__(self) -> Generator[str, None, None]:
        """Iterate over the keys of the Frame."""
        yield from self.raw.values.keys()
        yield from self.raw.arrays.keys()

    def __getitem__(self, k: Union[str, FrameKey]) -> Any:
        if isinstance(k, FrameKey):
            return k.get(self)
        try:
            return get_frame_key(k).get(self)
        except KeyError as e:
            if k in self.raw.values:
                return value_to_object(self.raw.values[k])
            if k in self.raw.arrays:
                arr = self.raw.arrays[k]
                if self.raw.arrays[k].HasField("index_values"):
                    return np.array(arr.ListFields()[0][1].values, dtype=int)
                elif self.raw.arrays[k].HasField("float_values"):
                    return np.array(arr.ListFields()[0][1].values, dtype=float)
                elif self.raw.arrays[k].HasField("string_values"):
                    return np.array(arr.ListFields()[0][1].values, dtype=object)
            raise KeyError from e

    def __setitem__(self, key: Union[str, FrameKey], value: Any) -> None:
        if isinstance(key, FrameKey):
            key.set(self, value)
            return
        try:
            get_frame_key(key).set(self, value)
        except KeyError:
            self._set_from_type(key, value)

    def _set_from_type(self, key: str, value: Any) -> None:
        if isinstance(value, str):
            self.raw.values[key].string_value = str(value)
        elif isinstance(value, float):
            self.raw.values[key].number_value = float(value)
        elif isinstance(value, np.ndarray):
            if value.dtype == float:
                self.set_float_array(key, value)
            elif value.dtype == int:
                if all(i >= 0 for i in value):
                    self.set_index_array(key, value)
                else:
                    self.set_float_array(key, value)
            elif value.dtype == object:
                self.set_string_array(key, value)
            else:
                raise TypeError(f"Did not know how to serialize {value}.")
        else:
            raise TypeError(f"Did not know how to serialize {value}.")

    def __contains__(self, key: Any) -> bool:
        if isinstance(key, FrameKey):
            return key.key in self.arrays or key.key in self.values
        return key in self.arrays or key in self.values

    @property
    def positions(self) -> Vector3Array:
        return self.particle_positions  # type: ignore

    @property
    def masses(self) -> ScalarArray:
        return self.particle_masses

    @property
    def velocities(self) -> Vector3Array:  # type: ignore
        return self.particle_velocities

    @classmethod
    def deserialize(cls, value: Serializable) -> FrameData:
        frame = FrameData()
        try:
            for key, v in value.items():  # type: ignore
                frame[key] = v
        except AttributeError as e:
            raise TypeError(f"{value} cannot be deserialized into a FrameData.") from e
        return frame

    def serialize(self) -> Serializable:
        json: Dict[str, Serializable] = {}
        for key, value in self.items():
            json[key] = value
        return json


def _print_value(value: Any) -> str:
    if isinstance(value, np.ndarray) and len(value) > 3:
        with np.printoptions(precision=3, suppress=True):
            return f"[{value[0]}, {value[1]}, ... ({len(value)} item(s)]"
    return repr(value)
