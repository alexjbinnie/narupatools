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
from typing import Iterator, List, Mapping, Optional

import numpy as np
import tables
from MDAnalysis.topology.tables import SYMB2Z
from infinite_sets import InfiniteSet
from narupa.trajectory import FrameData
from tables import File, Group, NoSuchNodeError

from narupatools.frame import (
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
from narupatools.frame.frame_source import FrameSource, TrajectorySource
from narupatools.physics.typing import IntArray, ScalarArray, Vector3Array


class InteractionView:
    """View of a specific interaction of a HDF5 trajectory."""

    def __init__(self, source: HDF5Trajectory, interaction: Group):
        self._source = source
        self._interaction = interaction
        self._start_index = interaction._v_attrs["startIndex"]
        try:
            self._end_index = int(interaction._v_attrs["endIndex"])
        except KeyError:
            self._end_index = int(interaction.frameIndex[-1])
        self._type = interaction._v_attrs["type"]
        self._indices = np.array(interaction.indices, dtype=int)

    @property
    def type(self) -> str:
        """Type of the interaction."""
        return self._type  # type: ignore

    @property
    def start_frame_index(self) -> int:
        """Index of the start frame."""
        return self._start_index  # type: ignore

    @property
    def end_frame_index(self) -> int:
        """Index of the end frame."""
        return self._end_index

    @property
    def start_time(self) -> float:
        """Start time of the interaction in picoseconds."""
        return self._source._times[self._start_index]  # type: ignore

    @property
    def end_time(self) -> float:
        """End time of the interaction in picoseconds."""
        return self._source._times[self._end_index]  # type: ignore

    @property
    def duration(self) -> float:
        """Duration of the interaction in picoseconds."""
        return self.end_time - self.start_time

    @property
    def indices(self) -> np.ndarray:
        """Indices of the particles affected by the interaction."""
        return self._indices

    @property
    def forces(self) -> Vector3Array:
        """
        Forces on the particles affected by the interaction.

        The forces are in kilojoules per mole per nanometer.
        """
        return np.array(self._interaction.forces)

    @property
    def potential_energies(self) -> ScalarArray:
        """
        Potential energies of the interaction at each frame.

        The energies are in kilojoules per mole.
        """
        return np.array(self._interaction.potentialEnergy)

    @property
    def interaction_positions(self) -> ScalarArray:
        """Position of the interaction at each frame affected by it, in nanometers."""
        return np.array(self._interaction.position)

    @property
    def interaction_scales(self) -> ScalarArray:
        """Scale of the interaction at each frame affected by it."""
        return np.array(self._interaction.scale)

    @property
    def frame_indices(self) -> IntArray:
        """Indices of the frames affected by the interaction."""
        return np.array(self._interaction.frameIndex, dtype=int)

    @property
    def frame_range(self) -> range:
        """Range of frame indices covered by this interaction."""
        return range(self.frame_indices[0], self.frame_indices[-1] + 1)

    def calculate_work(self) -> float:
        r"""
        Calculate the work done by an interaction, in kilojoules per mole.

        The work :math:`W` is given by:

        .. math:: W = \sum_i \int_{t_0}^{t_1} F_i(t) \cdot d x_i(t)

        where the sum over :math:`i` is over the particles affected by the interaction.
        In the case of a discretized trajectory, the trapezoidal rule for evaluating
        integrals is used:

        .. math:: W = \sum_i \sum_{n=n_0}^{n_1} \frac{1}{2} (F_i(t + n \delta t) +
                  F_i(t + (n+1) \delta t) \cdot (x_i(t + (n+1) \delta t) -
                  x_i (t + n \delta t))

        where now the second sum is over adjacent timesteps.
        """
        work = 0.0
        for rel_frame, abs_frame in enumerate(
            range(self._start_index, self._end_index)
        ):
            f0 = self._interaction.forces[rel_frame]
            f1 = self._interaction.forces[rel_frame + 1]
            s0 = self._source._positions[abs_frame]
            s1 = self._source._positions[abs_frame + 1]
            work_this_step = 0.0
            for rel_index, abs_index in enumerate(self.indices):
                # Use trapezoidal rule to calculate single step of integral F.dS
                F = 0.5 * (f0[rel_index] + f1[rel_index])
                dS = s1[abs_index] - s0[abs_index]
                work_this_step += np.dot(F, dS)
            work += work_this_step
        return work

    def calculate_cumulative_work(self) -> ScalarArray:
        """
        Calculate the cumulative work applied by this interaction for each frame.

        The size of this array is the same size as the number of frames this interaction
        was involved in.

        The work returned is in kilojoules per mole.
        """
        work = np.zeros(len(self.frame_indices))
        for rel_frame, abs_frame in enumerate(
            range(self._start_index, self._end_index)
        ):
            f0 = self._interaction.forces[rel_frame]
            f1 = self._interaction.forces[rel_frame + 1]
            s0 = self._source._positions[abs_frame]
            s1 = self._source._positions[abs_frame + 1]
            work_this_step = 0.0
            for rel_index, abs_index in enumerate(self.indices):
                # Use trapezoidal rule to calculate single step of integral F.dS
                F = 0.5 * (f0[rel_index] + f1[rel_index])
                dS = s1[abs_index] - s0[abs_index]
                work_this_step += np.dot(F, dS)
            work[rel_frame + 1] = work[rel_frame] + work_this_step
        return work

    def calculate_power(self) -> float:
        """
        Calculate the average power put in by this interaction.

        The power is in kilojoules per mole per picosecond.
        """
        work = self.calculate_work()
        return work / self.duration

    def __repr__(self) -> str:
        return "<InteractionView>"


class InteractionsView(Mapping[str, InteractionView]):
    """View of the interactions of a HDF5 trajectory."""

    def __init__(self, source: HDF5Trajectory, interactions: Optional[Group] = None):
        self._source = source
        self._interactions = interactions

    def __len__(self) -> int:
        if self._interactions is None:
            return 0
        return self._interactions._v_nchildren

    def __getitem__(self, item: str) -> InteractionView:
        if self._interactions is None:
            raise KeyError
        return InteractionView(self._source, self._interactions._v_children[item])

    def __iter__(self) -> Iterator[str]:
        if self._interactions is None:
            yield from ()
            return
        for child in self._interactions._f_iter_nodes():
            yield child._v_name

    def calculate_work(self) -> float:
        """Calculate the work done by all interactions, in kilojoules per mole."""
        if self._interactions is None:
            return 0.0
        return sum(interaction.calculate_work() for interaction in self.values())

    def calculate_cumulative_work(self) -> ScalarArray:
        """Calculate the work done by all interactions, in kilojoules per mole."""
        accum_work = np.zeros(len(self._source))
        if self._interactions is None:
            return accum_work
        for interaction in self.values():
            interaction_frames = interaction.frame_indices
            interaction_works = interaction.calculate_cumulative_work()
            accum_work[interaction_frames] += interaction_works
            accum_work[interaction_frames[-1] + 1 :] += interaction_works[-1]
        return accum_work

    @property
    def potential_energies(self) -> ScalarArray:
        """
        Potential energies for each frame due to interactions.

        These energies are in kilojoules per mole.
        """
        potential_energies = np.zeros(len(self._source))
        for interaction in self.values():
            potential_energies[
                interaction.frame_indices
            ] += interaction.potential_energies
        return potential_energies

    @property
    def forces(self) -> Vector3Array:
        """
        Forces on each particle for each frame due to all interactions.

        These forces are in kilojoules per mole per nanometer.
        """
        forces = np.zeros((len(self._source), 3))
        for interaction in self.values():
            forces[interaction.frame_indices] += interaction.forces
        return forces

    def __repr__(self) -> str:
        return f"<InteractionsView {len(self)} interaction(s)>"


class HDF5Topology(FrameSource):
    """Topology read from a HDF5 trajectory."""

    def __init__(self) -> None:
        self.atoms: List[HDF5Topology.Atom] = []
        self.residues: List[HDF5Topology.Residue] = []
        self.chains: List[HDF5Topology.Chain] = []
        self.bonds: List[List[int]] = []

    def _add_chain(self) -> HDF5Topology.Chain:
        chain = HDF5Topology.Chain(len(self.chains))
        self.chains.append(chain)
        return chain

    def _add_residue(self, chain: HDF5Topology.Chain) -> HDF5Topology.Residue:
        residue = HDF5Topology.Residue(len(self.residues), chain)
        self.residues.append(residue)
        return residue

    def _add_atom(self, residue: HDF5Topology.Residue) -> HDF5Topology.Atom:
        atom = HDF5Topology.Atom(len(self.atoms), residue)
        self.atoms.append(atom)
        return atom

    class Atom:
        """Atom in a HDF5 topology."""

        def __init__(self, index: int, residue: HDF5Topology.Residue):
            self.index = index
            self.residue = residue
            self.name: Optional[str] = None
            self.element: Optional[str] = None

        @property
        def atomic_number(self) -> Optional[int]:
            """Atomic number of the atom."""
            if self.element is not None:
                return SYMB2Z[self.element]
            return None

    class Residue:
        """Residue in a HDF5 topology."""

        def __init__(self, index: int, chain: HDF5Topology.Chain):
            self.index = index
            self.chain = chain
            self.name: Optional[str] = None

    class Chain:
        """Chain in a HDF5 topology."""

        def __init__(self, index: int):
            self.index = index
            self.name: Optional[str] = None

    @classmethod
    def from_string(cls, string: str) -> HDF5Topology:
        """Create a HDF5 topology from a JSON string."""
        content = json.loads(string)
        topology = cls()
        if "bonds" in content:
            topology.bonds = content["bonds"]
        if "chains" in content:
            for chain_json in content["chains"]:
                chain = topology._add_chain()
                if "name" in chain_json:
                    chain.name = chain_json["name"]
                if "residues" in chain_json:
                    for residue_json in chain_json["residues"]:
                        residue = topology._add_residue(chain)
                        if "name" in residue_json:
                            residue.name = residue_json["name"]
                        if "atoms" in residue_json:
                            for atom_json in residue_json["atoms"]:
                                atom = topology._add_atom(residue)
                                if "element" in atom_json:
                                    atom.element = atom_json["element"]
                                if "name" in atom_json:
                                    atom.name = atom_json["name"]
        return topology

    def get_frame(self, *, fields: InfiniteSet[str]) -> FrameData:  # noqa: D102
        frame = FrameData()
        if ParticleNames.key in fields:
            ParticleNames.set(frame, [atom.name or "" for atom in self.atoms])
        if ParticleElements.key in fields:
            ParticleElements.set(
                frame, [atom.atomic_number or 0 for atom in self.atoms]
            )
        if ParticleResidues.key in fields:
            ParticleResidues.set(frame, [atom.residue.index for atom in self.atoms])
        if ParticleCount.key in fields:
            ParticleCount.set(frame, len(self.residues))
        if ResidueNames.key in fields:
            ResidueNames.set(frame, [residue.name or "" for residue in self.residues])
        if ResidueChains.key in fields:
            ResidueChains.set(frame, [residue.chain.index for residue in self.residues])
        if ResidueCount.key in fields:
            ResidueCount.set(frame, len(self.residues))
        if ChainNames.key in fields:
            ChainNames.set(frame, [chain.name or "" for chain in self.chains])
        if ChainCount.key in fields:
            ChainCount.set(frame, len(self.chains))
        if BondCount.key in fields:
            BondCount.set(frame, len(self.bonds))
        if BondPairs.key in fields:
            BondPairs.set(frame, self.bonds)
        return frame

    def __repr__(self) -> str:
        return (
            f"<HDF5Topology {len(self.atoms)} atom(s), {len(self.residues)} "
            f"residues(s), {len(self.chains)} chain(s), {len(self.bonds)} bond(s)>"
        )


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
        self._topology = HDF5Topology.from_string(self._file.root.topology[0])
        self._positions = self._file.root.coordinates
        with contextlib.suppress(NoSuchNodeError):
            self._velocities = self._file.root.velocities
        with contextlib.suppress(NoSuchNodeError):
            self._forces = self._file.root.forces
        with contextlib.suppress(NoSuchNodeError):
            self._times = self._file.root.time
        with contextlib.suppress(NoSuchNodeError):
            self._kinetic_energies = self._file.root.kineticEnergy
        with contextlib.suppress(NoSuchNodeError):
            self._potential_energies = self._file.root.potentialEnergy

        self._interactions: InteractionsView = InteractionsView(self)
        with contextlib.suppress(NoSuchNodeError):
            self._interactions = InteractionsView(self, self._file.root.interactions)

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
        return len(self._positions)

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
        frame = self._topology.get_frame(fields=fields)
        if ParticleCount.key in fields:
            ParticleCount.set(frame, len(self._positions[0]))
        if ParticlePositions.key in fields:
            ParticlePositions.set(frame, self._positions[index])
        if ParticleVelocities.key in fields:
            ParticleVelocities.set(frame, self._velocities[index])
        if ParticleForces.key in fields:
            ParticleForces.set(frame, self._forces[index])
        if PotentialEnergy.key in fields:
            PotentialEnergy.set(frame, self._potential_energies[index])
        if KineticEnergy.key in fields:
            KineticEnergy.set(frame, self._kinetic_energies[index])
        return frame

    def __repr__(self) -> str:
        return (
            f"<HDF5Trajectory {len(self)} frame(s), "
            f"{self._times[-1] - self._times[0]:.3f} ps long, {len(self.interactions)} "
            f"interactions(s), {len(self._topology.atoms)} atom(s), "
            f"{len(self._topology.residues)} residues(s), {len(self._topology.chains)} "
            f"chain(s), {len(self._topology.bonds)} bond(s)>"
        )
