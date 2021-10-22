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

"""Reader for the NarupaTools HDF5 format."""

from __future__ import annotations

import contextlib
import json
import re
from dataclasses import dataclass
from typing import Iterator, List, Mapping, Optional

import numpy as np
import tables
from infinite_sets import InfiniteSet
from MDAnalysis.topology.tables import SYMB2Z
from narupa.trajectory import FrameData
from tables import File, Group, NoSuchNodeError

from narupatools.frame._frame_source import FrameSource, TrajectorySource
from narupatools.frame.fields import (
    BondCount,
    BondPairs,
    ChainCount,
    ChainNames,
    KineticEnergy,
    ParticleCount,
    ParticleElements,
    ParticleForces,
    ParticleNames,
    ParticlePositions,
    ParticleResidues,
    ParticleVelocities,
    PotentialEnergy,
    ResidueChains,
    ResidueCount,
    ResidueNames,
)
from narupatools.physics.typing import IntArray, ScalarArray, Vector3Array









class HDF5Trajectory(TrajectorySource):
    """Reader for the NarupaTools HDF5 format."""

    def __init__(self, file: File):
        self._file = file
        conventions = re.split(",| ", self._file.root._v_attrs["conventions"])
        if "Pande" not in conventions:
            raise ValueError(
                "Can't read HDF5 trajectory - Pande not specified in conventions."
            )
        mdtraj_conv_version = self._file.root._v_attrs["conventionVersion"]
        if mdtraj_conv_version != "1.1":
            raise ValueError(
                f"Can't read HDF5 trajectory - MDTraj conventions {mdtraj_conv_version}"
                f" unsupported."
            )
        narupatools_conv_version = self._file.root._v_attrs[
            "narupatoolsConventionVersion"
        ]
        if narupatools_conv_version != "1.0":
            raise ValueError(
                f"Can't read HDF5 trajectory - NarupaTools conventions "
                f"{narupatools_conv_version} unsupported."
            )
        self._positions = self._file.root.coordinates
        self._n_frames, self._n_atoms, dim = self._positions.shape
        if dim != 3:
            raise ValueError("Trajectory does not contain 3D positions.")
        with contextlib.suppress(NoSuchNodeError):
            self._topology = HDF5Topology.from_string(self._file.root.topology[0])
            if len(self._topology.atoms) != self._n_atoms:
                raise ValueError("Topology has wrong number of atoms.")
        with contextlib.suppress(NoSuchNodeError):
            self._velocities = self._file.root.velocities
            if self._velocities.shape != self._positions.shape:
                raise ValueError("Position and velocity arrays have mismatched shapes.")
        with contextlib.suppress(NoSuchNodeError):
            self._forces = self._file.root.forces
            if self._forces.shape != self._positions.shape:
                raise ValueError("Position and forces arrays have mismatched shapes.")
        with contextlib.suppress(NoSuchNodeError):
            self._times = self._file.root.time
            if self._times.shape != (self._n_frames,):
                raise ValueError("Time array has mismatched shape.")
        with contextlib.suppress(NoSuchNodeError):
            self._kinetic_energies = self._file.root.kineticEnergy
            if self._kinetic_energies.shape != (self._n_frames,):
                raise ValueError("Kinetic energy array has mismatched shape.")
        with contextlib.suppress(NoSuchNodeError):
            self._potential_energies = self._file.root.potentialEnergy
            if self._potential_energies.shape != (self._n_frames,):
                raise ValueError("Potential energy array has mismatched shape.")

        self._interactions: InteractionsView = InteractionsView(self)
        with contextlib.suppress(NoSuchNodeError):
            self._interactions = InteractionsView(self, self._file.root.interactions)

    def close(self) -> None:
        """Close the trajectory and the underlying file."""
        self._file.close()

    @classmethod
    def load_file(cls, filename: str) -> HDF5Trajectory:
        """
        Load a trajectory that conforms to the narupatools HDF5 format.

        :param filename: Filename to load.
        """
        file = tables.open_file(filename, mode="r")
        return cls(file)

    @property
    def interactions(self) -> InteractionsView:
        """Interactions applied during this trajectory."""
        return self._interactions

    def __len__(self) -> int:
        return self._n_frames  # type: ignore

    @property
    def topology(self) -> HDF5Topology:
        """Topology stored with the trajectory."""
        return self._topology

    @property
    def positions(self) -> Vector3Array:
        """Positions of the particles for each frame, in nanometers."""
        return np.array(self._positions)

    @property
    def velocities(self) -> Vector3Array:
        """Velocities of the particles for each frame, in nanometers per picoseconds."""
        return np.array(self._velocities)

    @property
    def forces(self) -> Vector3Array:
        """Forces on the particles for each frame, in nanometers per picoseconds."""
        return np.array(self._forces)

    @property
    def times(self) -> ScalarArray:
        """Time of each frame, in picoseconds."""
        return np.array(self._times)

    @property
    def kinetic_energies(self) -> ScalarArray:
        """Kinetic energy of each frame, in kilojoules per mole."""
        return np.array(self._kinetic_energies)

    @property
    def potential_energies(self) -> ScalarArray:
        """Potential energy of each frame, in kilojoules per mole."""
        return np.array(self._potential_energies)

    def get_frame(  # noqa: D102
        self, *, index: int, fields: InfiniteSet[str]
    ) -> FrameData:
        frame = FrameData()
        with contextlib.suppress(AttributeError):
            frame = self._topology.get_frame(fields=fields)
        if ParticleCount.key in fields:
            frame[ParticleCount] = len(self._positions[0])
        if ParticlePositions.key in fields:
            frame[ParticlePositions] = self._positions[index]
        if ParticleVelocities.key in fields:
            with contextlib.suppress(AttributeError):
                frame[ParticleVelocities] = self._velocities[index]
        if ParticleForces.key in fields:
            with contextlib.suppress(AttributeError):
                frame[ParticleForces] = self._forces[index]
        if PotentialEnergy.key in fields:
            with contextlib.suppress(AttributeError):
                frame[PotentialEnergy] = self._potential_energies[index]
        if KineticEnergy.key in fields:
            with contextlib.suppress(AttributeError):
                frame[KineticEnergy] = self._kinetic_energies[index]
        return frame

    def __repr__(self) -> str:
        return (
            f"<HDF5Trajectory {len(self)} frame(s), "
            f"{self._times[-1] - self._times[0]:.3f} ps long, {len(self.interactions)} "
            f"interactions(s), {len(self._topology.atoms)} atom(s), "
            f"{len(self._topology.residues)} residues(s), {len(self._topology.chains)} "
            f"chain(s), {len(self._topology.bonds)} bond(s)>"
        )
