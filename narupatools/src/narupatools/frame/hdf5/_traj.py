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

from __future__ import annotations

import contextlib
import os
import warnings
from tempfile import NamedTemporaryFile
from time import monotonic
from typing import Any, Generator, Iterator, Mapping, Optional

import numpy as np
import numpy.typing as npt
import tables
from infinite_sets import InfiniteSet, everything
from narupa.trajectory import FrameData
from tables import File, Group, NoSuchNodeError

import narupatools
from narupatools.frame import (
    DynamicStructureMethods,
    KineticEnergy,
    ParticleForces,
    ParticlePositions,
    ParticleVelocities,
    PotentialEnergy,
    StateData,
    TrajectorySource,
)
from narupatools.imd import Interaction, InteractiveSimulationDynamics
from narupatools.physics.energy import cumulative_work, total_work
from narupatools.physics.typing import ScalarArray, Vector3Array

from ._descriptor import HDF5AppendableArray, HDF5Attribute
from ._object import _HDF5EditableObject
from ._topology import HDF5Topology
from ._utils import generate_topology


class HDF5InteractionParameters(_HDF5EditableObject):
    """Interaction parameters stored in a trajectory."""

    def __init__(self, interaction_view: HDF5Interaction):
        self._interaction = interaction_view

    @property
    def writable(self) -> bool:  # noqa: D102
        return self._interaction.writable

    @property
    def hdf5_group(self) -> Group:  # noqa: D102
        return self._interaction.hdf5_group

    @property
    def n_atoms(self) -> int:
        """Number of atoms involved in the interaction."""
        return self._interaction.n_atoms

    _position = HDF5AppendableArray(
        h5_name="position", title="Interaction Position", shape=(3,), units="nanometers"
    )
    position = _position.as_numpy()

    _scale = HDF5AppendableArray(
        h5_name="scale", title="Interaction Scale", shape=(), units=None
    )
    scale = _scale.as_numpy()

    _force = HDF5AppendableArray(
        h5_name="force",
        title="Interaction Force",
        shape=(3,),
        units="kilojoules/mole/angstrom",
    )
    force = _force.as_numpy()

    def save_interaction(self, *, interaction: Interaction) -> None:
        """Save interaction parameters to file."""
        if hasattr(interaction, "position"):
            self._position.append([interaction.position])  # type: ignore[attr-defined]
        if hasattr(interaction, "scale"):
            self._scale.append([interaction.scale])  # type: ignore[attr-defined]
        if hasattr(interaction, "force"):
            self._force.append([interaction.force])  # type: ignore[attr-defined]


class HDF5Interaction(_HDF5EditableObject):
    """View of a specific interaction of a HDF5 trajectory."""

    start_frame_index = HDF5Attribute("startIndex")
    end_frame_index = HDF5Attribute("endIndex")
    interaction_type = HDF5Attribute("type")
    """Type of the interaction."""

    def __init__(self, trajectory: HDF5Trajectory, interaction: Group):
        self._trajectory = trajectory
        """Trajectory the interaction belongs to."""
        self._group = interaction
        """HDF5 group of this interaction."""
        self._particle_indices: Optional[np.ndarray] = None
        self.parameters = HDF5InteractionParameters(self)

    @property
    def n_atoms(self) -> int:
        """Number of atoms in the interaction."""
        return len(self.particle_indices)

    @property
    def writable(self) -> bool:  # noqa: D102
        return self._trajectory.writable

    @property
    def hdf5_group(self) -> Group:  # noqa: D102
        return self._group

    @classmethod
    def create(
        cls,
        *,
        key: str,
        trajectory: HDF5Trajectory,
        interaction_group: Group,
        interaction: Interaction,
        frame_index: int,
    ) -> HDF5Interaction:
        """Create a new interaction."""
        group = trajectory._file.create_group(
            where=interaction_group, name=key, title="Interaction"
        )

        view = HDF5Interaction(trajectory, group)

        view.interaction_type = interaction.interaction_type
        view.start_frame_index = frame_index
        view._set_particle_indices(interaction.particle_indices)

        return view

    @property
    def _end_index(self) -> int:
        try:
            return int(self.end_frame_index)
        except KeyError:
            return int(self.frame_indices[-1])

    _forces = HDF5AppendableArray(
        h5_name="forces",
        title="Interactive Forces",
        per_atom=True,
        shape=(3,),
        units="kilojoules/(mole * nanometer",
    )
    forces = _forces.as_numpy()

    _potential_energies = HDF5AppendableArray(
        h5_name="potentialEnergy",
        title="Interactive Potential Energy",
        shape=(),
        units="kilojoules/mole",
    )
    potential_energies = _potential_energies.as_numpy()

    _frame_indices = HDF5AppendableArray(
        h5_name="frameIndex", title="Frame Index", shape=(), units=None
    )
    frame_indices = _frame_indices.as_numpy(dtype=int)

    def _set_particle_indices(self, particle_indices: npt.ArrayLike) -> None:
        if not self.writable:
            raise ValueError("Trajectory is not writable")

        self._particle_indices = np.asarray(particle_indices, dtype=int)

        self._trajectory._file.create_array(
            where=self._group,
            name="indices",
            title="Particle Indices",
            obj=self._particle_indices,
        )

    @property
    def particle_indices(self) -> npt.NDArray[np.int_]:
        """Particle indices involved in the interaction."""
        if self._particle_indices is None:
            if hasattr(self._group, "indices"):
                self._particle_indices = np.array(self._group.indices)
            else:
                raise AttributeError
        return self._particle_indices

    @property
    def indices(self) -> npt.NDArray[np.int_]:
        """
        Particle indices involved in the interaction.

        This is deprecated, and particle_indices should be used instead.
        """
        warnings.simplefilter("always", DeprecationWarning)
        warnings.warn(
            f"Call to deprecated function {HDF5Interaction.indices.__name__}.",
            category=DeprecationWarning,
            stacklevel=2,
        )
        warnings.simplefilter("default", DeprecationWarning)
        return self.particle_indices

    def save_interaction(
        self,
        *,
        interaction: Interaction,
        frame_index: int,
    ) -> None:
        """
        Save an IMD interaction to the trajectory.

        :param frame_index: Index of the frame in the main trajectory this interaction
                            corresponds to.
        :param interaction: Interaction to save.
        """
        self._potential_energies.append([interaction.potential_energy])
        self._forces.append([interaction.forces])
        self.parameters.save_interaction(interaction=interaction)
        self._frame_indices.append([frame_index])

    @property
    def start_time(self) -> float:
        """Start time of the interaction in picoseconds."""
        return self._trajectory._times[self.start_frame_index]  # type: ignore

    @property
    def end_time(self) -> float:
        """End time of the interaction in picoseconds."""
        return self._trajectory._times[self.end_frame_index]  # type: ignore

    @property
    def duration(self) -> float:
        """Duration of the interaction in picoseconds."""
        return self.end_time - self.start_time

    @property
    def frame_range(self) -> range:
        """Range of frame indices covered by this interaction."""
        return range(self.frame_indices[0], self.frame_indices[-1] + 1)

    @property
    def frame_slice(self) -> slice:
        """Slice of frame indices covered by this interaction."""
        return slice(self.frame_indices[0], self.frame_indices[-1] + 1)

    def calculate_per_particle_work(self) -> np.ndarray:
        """Calculate the per-particle work done by all interactions."""
        return total_work(  # type: ignore
            forces=self.forces,
            positions=self._trajectory.positions[
                ..., self.frame_slice, self.particle_indices, :
            ],
            time_axis=-3,
        )

    def calculate_cumulative_work(self) -> np.ndarray:
        """Calculate the cumulative work done by the interaction at each timestep."""
        return cumulative_work(  # type: ignore
            forces=self.forces,
            positions=self._trajectory.positions[
                ..., self.frame_slice, self.particle_indices, :
            ],
            time_axis=-3,
        ).sum(axis=-1)

    def calculate_work(self) -> np.ndarray:
        """Calculate the cumulative work done by the interaction at each timestep."""
        return total_work(  # type: ignore
            forces=self.forces,
            positions=self._trajectory.positions[
                ..., self.frame_slice, self.particle_indices, :
            ],
            time_axis=-3,
        ).sum(axis=-1)

    def __str__(self) -> str:
        return f"<InteractionView indices={self.particle_indices} start_index={self.start_frame_index} end_index={self.end_frame_index} type={self.interaction_type}>"


class InteractionsView(Mapping[str, HDF5Interaction]):
    """View of the interactions of a HDF5 trajectory."""

    def __init__(self, trajectory: HDF5Trajectory):
        self._trajectory = trajectory
        self._interactions_group: Optional[Group] = None

    @property
    def writable(self) -> bool:  # noqa: D102
        return self._trajectory.writable

    @property
    def _interactions(self) -> Group:
        if self._interactions_group is None:
            if hasattr(self._trajectory._file.root, "interactions"):
                self._interactions_group = self._trajectory._file.root.interactions
            elif self.writable:
                self._interactions_group = self._trajectory._file.create_group(
                    where=self._trajectory._file.root,
                    name="interactions",
                    title="Interactions",
                )
            else:
                raise AttributeError
        return self._interactions_group

    def __len__(self) -> int:
        try:
            return self._interactions._v_nchildren
        except AttributeError:
            return 0

    def __getitem__(self, item: str) -> HDF5Interaction:
        try:
            interactions = self._interactions
        except AttributeError as e:
            raise KeyError from e
        try:
            return HDF5Interaction(self._trajectory, interactions._f_get_child(item))
        except NoSuchNodeError as e:
            raise KeyError from e

    def __iter__(self) -> Iterator[str]:
        try:
            for child in self._interactions._f_iter_nodes():
                yield child._v_name
        except AttributeError:
            yield from ()
            return

    def save_interaction(
        self,
        *,
        key: str,
        interaction: Interaction,
        frame_index: int,
    ) -> None:
        """
        Save an IMD interaction to the trajectory.

        :param frame_index: Index of the frame in the main trajectory this interaction
                            corresponds to.
        :param key: Unique key of interaction.
        :param interaction: Interaction to save.
        """
        if not self.writable:
            raise ValueError("Not writable")
        if key not in self:
            HDF5Interaction.create(
                key=key,
                trajectory=self._trajectory,
                interaction_group=self._interactions,
                interaction=interaction,
                frame_index=frame_index,
            )
        self[key].save_interaction(interaction=interaction, frame_index=frame_index)

    def end_interaction(self, *, key: str, frame_index: int) -> None:
        """
        Mark an interaction as having finished.

        :param key: Key of the interaction.
        :param frame_index: Current frame of the simulation.
        """
        if key not in self:
            return
        self[key].end_frame_index = frame_index

    @property
    def potential_energies(self) -> ScalarArray:
        """
        Potential energies for each frame due to interactions.

        These energies are in kilojoules per mole.
        """
        potential_energies = np.zeros(self._trajectory.n_frames)
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
        forces = np.zeros((self._trajectory.n_frames, self._trajectory.n_atoms, 3))
        for interaction in self.values():
            forces[
                interaction.frame_slice, interaction.particle_indices
            ] += interaction.forces
        return forces

    def calculate_per_particle_work(self) -> np.ndarray:
        """Calculate the per-particle work done by all interactions."""
        return total_work(  # type: ignore
            forces=self.forces, positions=self._trajectory.positions, time_axis=-3
        )

    def calculate_cumulative_work(self) -> np.ndarray:
        """Calculate the cumulative work done by all interactions at each timestep."""
        return cumulative_work(  # type: ignore
            forces=self.forces, positions=self._trajectory.positions, time_axis=-3
        ).sum(axis=-1)

    def calculate_work(self) -> np.ndarray:
        """Calculate the work done by all interactions at each timestep."""
        return total_work(  # type: ignore
            forces=self.forces, positions=self._trajectory.positions, time_axis=-3
        ).sum(axis=-1)

    def __repr__(self) -> str:
        return f"<InteractionsView {len(self)} interaction(s)>"


class HDF5Trajectory(DynamicStructureMethods, TrajectorySource, _HDF5EditableObject):
    """Trajectory stored as an HDF5 trajectory."""

    application = HDF5Attribute("application")
    convention_version = HDF5Attribute("conventionVersion")
    narupatools_convention_version = HDF5Attribute("narupatoolsConventionVersion")
    program = HDF5Attribute("program")
    conventions = HDF5Attribute("conventions")
    program_version = HDF5Attribute("programVersion")
    author = HDF5Attribute("author")

    _positions = HDF5AppendableArray(
        h5_name="coordinates",
        title="Atomic Coordinates",
        shape=(3,),
        per_atom=True,
        units="nanometers",
    )
    positions = _positions.as_numpy()

    _velocities = HDF5AppendableArray(
        h5_name="velocities",
        title="Atomic Velocities",
        shape=(3,),
        per_atom=True,
        units="nanometer/picosecond",
    )
    velocities = _velocities.as_numpy()

    _forces = HDF5AppendableArray(
        h5_name="forces",
        title="Atomic Forces",
        shape=(3,),
        per_atom=True,
        units="nanometers",
    )
    forces = _forces.as_numpy()

    _kinetic_energies = HDF5AppendableArray(
        h5_name="kineticEnergy",
        title="Kinetic Energy",
        shape=(),
        units="kilojoules/mole",
    )
    kinetic_energies = _kinetic_energies.as_numpy()

    _potential_energies = HDF5AppendableArray(
        h5_name="potentialEnergy",
        title="Potential Energy",
        shape=(),
        units="kilojoules/mole",
    )
    potential_energies = _potential_energies.as_numpy()

    _times = HDF5AppendableArray(
        h5_name="time", title="Time", shape=(), units="picoseconds"
    )
    times = _times.as_numpy()

    _realtimes = HDF5AppendableArray(
        h5_name="realtime", title="Real Time", shape=(), units="seconds"
    )
    realtimes = _realtimes.as_numpy()

    def __len__(self) -> int:
        return self.n_frames

    @classmethod
    def load_file(cls, filename: str) -> HDF5Trajectory:
        """Load a HDF5 trajectory from a file."""
        file = tables.open_file(filename, mode="r")
        obj = cls(file, False)

        if "Pande" not in obj.conventions:
            raise ValueError(
                "Can't read HDF5 trajectory - Pande not specified in conventions."
            )
        if obj.convention_version != "1.1":
            raise ValueError(
                f"Can't read HDF5 trajectory - MDTraj conventions {obj.convention_version}"
                f" unsupported."
            )
        if obj.narupatools_convention_version != "1.0":
            raise ValueError(
                f"Can't read HDF5 trajectory - NarupaTools conventions "
                f"{obj.narupatools_convention_version} unsupported."
            )

        return obj

    @classmethod
    def new_file(
        cls,
        filename: str,
        *,
        title: Optional[str] = None,
        author: Optional[str] = None,
        overwrite_existing: bool = False,
        expected_frames: Optional[int] = None,
    ) -> HDF5Trajectory:
        """Create a trajectory streamed to a file."""
        if not overwrite_existing and os.path.exists(filename):
            raise FileExistsError(filename)
        file = tables.open_file(filename, mode="w", title=title or "")
        obj = cls(file, True)

        obj.author = author
        obj.application = "narupatools"
        obj.convention_version = "1.1"
        obj.narupatools_convention_version = "1.0"
        obj.program = "narupatools"
        obj.conventions = "Pande,NarupaTools"
        obj.program_version = narupatools.__version__
        obj.expected_frames = expected_frames

        return obj

    @classmethod
    def in_memory(
        cls,
        *,
        title: Optional[str] = None,
        author: Optional[str] = None,
        expected_frames: Optional[int] = None,
    ) -> HDF5Trajectory:
        """Create a trajectory stored in memory."""
        with NamedTemporaryFile() as tmp:
            filename = tmp.name + ".h5"
        file = tables.open_file(
            filename,
            mode="w",
            title=title or "",
            driver="H5FD_CORE",
            driver_core_backing_store=0,
        )

        obj = cls(file, True)

        obj.author = author
        obj.application = "narupatools"
        obj.convention_version = "1.1"
        obj.narupatools_convention_version = "1.0"
        obj.program = "narupatools"
        obj.conventions = "Pande,NarupaTools"
        obj.program_version = narupatools.__version__
        obj.expected_frames = expected_frames

        return obj

    def __init__(self, file: File, writable: bool):
        self._file = file
        self._writable = writable
        self._start_realtime: Optional[float] = None
        self._topology: Optional[HDF5Topology] = None
        self.expected_frames: Optional[int] = None
        self.interactions = InteractionsView(self)

    @property
    def writable(self) -> bool:  # noqa: D102
        return self._writable

    @property
    def hdf5_group(self) -> Group:  # noqa: D102
        return self._file.root

    @property
    def n_atoms(self) -> int:
        """Number of atoms in the trajectory."""
        if hasattr(self, "_n_atoms"):
            return self._n_atoms  # type: ignore
        self._n_atoms = int(self._positions.shape[1])
        return self._n_atoms

    @property
    def n_frames(self) -> int:
        """Number of frames in the trajectory."""
        return self.positions.shape[0]

    def save_frame(self, *, frame: FrameData, time: float) -> None:
        """
        Save a trajectory frame.

        :param time: Time in picoseconds.
        """
        if not self._writable:
            raise ValueError("Trajectory is read only")

        positions = frame[ParticlePositions]

        if not hasattr(self, "_n_atoms"):
            self._n_atoms = len(positions)

        self._positions.append(positions[np.newaxis])
        if ParticleVelocities in frame:
            self._velocities.append(frame[ParticleVelocities][np.newaxis])
        if ParticleForces in frame:
            self._forces.append(frame[ParticleForces][np.newaxis])
        self._times.append(np.array([time]))
        if self._start_realtime is None:
            self._start_realtime = monotonic()
            self._realtimes.append(np.array([0]))
        else:
            self._realtimes.append(np.array([monotonic() - self._start_realtime]))

        self._potential_energies.append(np.array([frame[PotentialEnergy]]))
        self._kinetic_energies.append(np.array([frame[KineticEnergy]]))

    def save_topology(self, frame: FrameData) -> None:
        """
        Save the topology of the given frame to the file.

        :param frame: Frame to get topology from.
        """
        if not self._writable:
            raise ValueError("Trajectory is read only")

        raw_string = generate_topology(frame)
        self._file.create_array(
            self._file.root, name="topology", title="Topology", obj=[raw_string]
        )
        self._topology = HDF5Topology.from_string(raw_string.decode("ascii"))

    def flush(self) -> None:
        """Flush outstanding changes to file."""
        self._file.flush()

    @property
    def topology(self) -> HDF5Topology:
        """Topology of the system."""
        if self._topology is None:
            if "topology" in self._file.root:
                self._topology = HDF5Topology.from_string(
                    self._file.root.topology[0].decode("ascii")
                )
            else:
                raise AttributeError
        return self._topology

    @property
    def masses(self) -> np.ndarray:
        """Masses of each atom in the system, in daltons."""
        return self.topology.masses

    def get_frame(  # noqa: D102
        self, *, index: int, fields: InfiniteSet[str] = everything()
    ) -> FrameData:
        frame = FrameData()
        with contextlib.suppress(AttributeError):
            frame = self.topology.get_frame(fields=fields)
        if ParticlePositions in fields:
            frame[ParticlePositions] = self._positions[index]
        if ParticleVelocities in fields:
            with contextlib.suppress(AttributeError):
                frame[ParticleVelocities] = self._velocities[index]
        if ParticleForces in fields:
            with contextlib.suppress(AttributeError):
                frame[ParticleForces] = self._forces[index]
        if PotentialEnergy in fields:
            with contextlib.suppress(AttributeError):
                frame[PotentialEnergy] = self._potential_energies[index]
        if KineticEnergy in fields:
            with contextlib.suppress(AttributeError):
                frame[KineticEnergy] = self._kinetic_energies[index]
        return frame

    def save_to_file(self, filename: str, overwrite_existing: bool = False) -> None:
        """Save the trajectory to the given filename."""
        self._file.flush()
        self._file.copy_file(filename, overwrite_existing)

    def close(self) -> None:
        """Flush any outstanding changes to file and close the file."""
        self._file.flush()
        self._file.close()

    @classmethod
    @contextlib.contextmanager
    def record(
        cls,
        dynamics: InteractiveSimulationDynamics,
        *,
        filename: Optional[str] = None,
        expected_frames: Optional[int] = None,
        flush_every: int = 1,
        title: Optional[str] = None,
        close_file_after: bool = False,
        overwrite_existing: bool = False,
    ) -> Generator[HDF5Trajectory, None, None]:
        """Record dynamics to a single trajectory, stopping if it is reset."""
        if filename is not None:
            traj = cls.new_file(
                filename=filename,
                expected_frames=expected_frames,
                title=title,
                overwrite_existing=overwrite_existing,
            )
        else:
            traj = cls.in_memory(expected_frames=expected_frames, title=title)

        has_logged_initial = False
        frame_index = 0

        def log_initial(**kwargs: Any) -> None:
            nonlocal traj
            nonlocal has_logged_initial

            if has_logged_initial:
                # Interactions that have appeared should be registered now.
                for interaction in dynamics.imd.current_interactions.values():
                    if interaction.key not in traj.interactions:
                        traj.interactions.save_interaction(
                            key=interaction.key,
                            interaction=interaction,
                            frame_index=frame_index - 1,
                        )
            else:
                log_step(**kwargs)
                traj.save_topology(dynamics.get_frame(fields=everything()))
                has_logged_initial = True

        def log_step(**kwargs: Any) -> None:
            nonlocal traj
            nonlocal frame_index

            fields = {
                ParticlePositions,
                ParticleVelocities,
                ParticleForces,
                KineticEnergy,
                PotentialEnergy,
            }

            state_data = StateData()

            frame = dynamics.get_frame(fields=fields, existing=state_data)  # type: ignore

            traj.save_frame(
                frame=frame,
                time=dynamics.elapsed_time,
            )
            for interaction in dynamics.imd.current_interactions.values():
                traj.interactions.save_interaction(
                    key=interaction.key,
                    interaction=interaction,
                    frame_index=frame_index,
                )

            frame_index += 1

            if frame_index % flush_every == 0:
                traj.flush()

        def log_end_interaction(key: str, **kwargs: Any) -> None:
            nonlocal traj
            nonlocal frame_index

            traj.interactions.end_interaction(key=key, frame_index=frame_index - 1)

        dynamics.on_pre_step.add_callback(log_initial, priority=-1000)
        dynamics.on_post_step.add_callback(log_step, priority=-1000)
        dynamics.imd.on_end_interaction.add_callback(log_end_interaction)

        try:
            yield traj
        finally:
            traj.flush()
            dynamics.on_pre_step.remove_callback(log_initial)
            dynamics.on_post_step.remove_callback(log_step)
            dynamics.imd.on_end_interaction.remove_callback(log_end_interaction)
            traj._writable = False
            if close_file_after:
                traj.close()

    def __str__(self) -> str:
        s = "<HDF5Trajectory"
        if not self._writable:
            s += " read-only"
        s += f" {len(self.positions)} frame(s)"
        s += f" {self.n_atoms} atom(s)"
        s += ">"
        return s
