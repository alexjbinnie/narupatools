from __future__ import annotations

import contextlib
import os
import re
from abc import abstractmethod
from tempfile import NamedTemporaryFile
from time import monotonic
from typing import Optional, Tuple, Mapping, Iterator, Protocol, Any

import numpy as np
import tables
from infinite_sets import InfiniteSet, everything
from narupa.trajectory import FrameData
from tables import File, EArray, Float32Atom, Group

import narupatools
from narupatools.core.dynamics import SimulationDynamics
from narupatools.frame import ParticlePositions, ParticleVelocities, ParticleForces, PotentialEnergy, \
    KineticEnergy, DynamicStructureMethods, TrajectorySource
from narupatools.frame.hdf5._topology import HDF5Topology
from narupatools.frame.hdf5._utils import generate_topology
from narupatools.imd import Interaction, InteractiveSimulationDynamics
from narupatools.physics.typing import Vector3Array, ScalarArray


class _HDF5EditableObject(Protocol):
    """Define an object backed by a HDF5 group."""

    @property
    @abstractmethod
    def hdf5_group(self) -> Group:
        pass

    @property
    @abstractmethod
    def writable(self) -> bool:
        pass

    def create_earray(
            self,
            *,
            name: str,
            title: str,
            shape: Tuple[int, ...],
            units: Optional[str] = None,
    ) -> EArray:
        array = self.hdf5_group._v_file.create_earray(
            self.hdf5_group,
            name=name,
            atom=Float32Atom(),
            shape=(0,) + shape,
            title=title,
            filters=tables.Filters(shuffle=True, complib="zlib", complevel=1),
        )
        if units is not None:
            array._v_attrs["units"] = units
        return array


class HDF5Attribute:

    def __init__(self, name: str):
        self._name = name

    def __get__(self, instance: _HDF5EditableObject, objtype):
        try:
            return instance.hdf5_group._v_attrs[self._name]
        except KeyError as e:
            raise AttributeError(f"Attribute {self._name} is not present.") from e

    def __set__(self, instance, value) -> None:
        if not instance.writable:
            raise ValueError("Trajectory is not writable.")
        instance.hdf5_group._v_attrs[self._name] = value




class HDF5AppendableArray:
    """
    Descriptor representing a HDF5 array that can be appended to.

    When first accessed, the following steps are taken:

    * If the array exists (due to being read from an existing file), it is returned.
    * If
    """
    def __init__(
            self,
            *,
            h5_name: str,
            title: str,
            shape: Tuple[int],
            per_atom: bool = False,
            units: str,
    ):
        self._h5_name = h5_name
        """Name of the array stored as an HDF5 EArray."""
        self._array_name = f"_{h5_name}_array"
        """Name to store the array as a python attribute."""
        self._title = title
        """Descriptive title of the array."""
        self._per_atom = per_atom
        """Should the size of the array be multiplied by the number of atoms?"""
        self._shape = shape
        """Shape of the numpy array."""
        self._units = units
        """Units of the array."""

    def __get__(self, instance: _HDF5EditableObject, objtype):
        try:
            return getattr(instance, self._array_name)
        except AttributeError:
            if self._h5_name in instance.hdf5_group:
                value = getattr(instance.hdf5_group, self._h5_name)
            elif instance.writable:
                shape = self._shape
                if self._per_atom:
                    shape = (instance.n_atoms, *shape)
                value = instance.create_earray(
                    name=self._h5_name,
                    title=self._title,
                    shape=shape,
                    units=self._units,
                )
            else:
                raise AttributeError
            setattr(instance, self._array_name, value)
            return value

    def as_numpy(self, dtype=float):
        def get(obj):
            return np.array(self.__get__(obj, type(obj)), dtype=dtype)

        return property(fget=get)


class HDF5InteractionParameters(_HDF5EditableObject):

    def __init__(self, interaction_view: HDF5Interaction):
        self._interaction = interaction_view

    @property
    def writable(self) -> bool:
        return self._interaction.writable

    @property
    def hdf5_group(self) -> Group:
        return self._interaction.hdf5_group

    _position = HDF5AppendableArray(h5_name="position", title="Interaction Position", shape=(3,), units="nanometers")
    position = _position.as_numpy()

    _scale = HDF5AppendableArray(h5_name="scale", title="Interaction Scale", shape=(), units=None)
    scale = _scale.as_numpy()

    _force = HDF5AppendableArray(h5_name="force", title="Interaction Force", shape=(3,), units="kilojoules/mole/angstrom")
    force = _force.as_numpy()

    def save_interaction(self, *, interaction: Interaction):
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
        self._particle_indices = None
        self.parameters = HDF5InteractionParameters(self)

    @property
    def n_atoms(self):
        return len(self.particle_indices)

    @property
    def writable(self):
        return self._trajectory.writable

    @property
    def hdf5_group(self) -> Group:
        return self._group

    @classmethod
    def create(cls, *, key: str, trajectory: HDF5Trajectory, interaction_group: Group, interaction: Interaction,
               frame_index: int):

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

    _forces = HDF5AppendableArray(h5_name="forces", title="Interactive Forces", per_atom=True, shape=(3,),
                                  units="kilojoules/(mole * nanometer")
    forces = _forces.as_numpy()

    _potential_energies = HDF5AppendableArray(h5_name="potentialEnergy", title="Interactive Potential Energy", shape=(),
                                              units="kilojoules/mole")
    potential_energies = _potential_energies.as_numpy()

    _frame_indices = HDF5AppendableArray(h5_name="frameIndex", title="Frame Index", shape=(), units=None)
    frame_indices = _frame_indices.as_numpy(dtype=int)

    def _set_particle_indices(self, particle_indices):
        if not self.writable:
            raise ValueError("Trajectory is not writable")
        self._trajectory._file.create_array(
            where=self._group,
            name="indices",
            title="Particle Indices",
            obj=list(particle_indices),
        )
        self._particle_indices = particle_indices

    @property
    def particle_indices(self):
        if self._particle_indices is None:
            if hasattr(self._group, "indices"):
                self._particle_indices = np.array(self._group.indices)
            else:
                raise AttributeError
        return self._particle_indices

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

    def __str__(self) -> str:
        return f"<InteractionView indices={self.particle_indices} start_index={self.start_frame_index} end_index={self.end_frame_index} type={self.interaction_type}>"


class InteractionsView(Mapping[str, HDF5Interaction]):
    """View of the interactions of a HDF5 trajectory."""

    def __init__(self, trajectory: HDF5Trajectory):
        self._trajectory = trajectory
        self._interactions_group = None

    @property
    def writable(self):
        return self._trajectory.writable

    @property
    def _interactions(self):
        if self._interactions_group is None:
            if hasattr(self._trajectory._file.root, 'interactions'):
                self._interactions_group = self._trajectory._file.root.interactions
            elif self.writable:
                self._interactions_group = self._trajectory._file.create_group(
                    where=self._trajectory._file.root, name="interactions", title="Interactions"
                )
            else:
                raise AttributeError
        return self._interactions_group

    def __len__(self) -> int:
        if self._interactions is None:
            return 0
        return self._interactions._v_nchildren

    def __getitem__(self, item: str) -> HDF5Interaction:
        if self._interactions is None:
            raise KeyError
        if hasattr(self._interactions, item):
            return HDF5Interaction(self._trajectory, getattr(self._interactions, item))
        else:
            raise KeyError

    def __iter__(self) -> Iterator[str]:
        if self._interactions is None:
            yield from ()
            return
        for child in self._interactions._f_iter_nodes():
            yield child._v_name

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
            HDF5Interaction.create(key=key, trajectory=self._trajectory, interaction_group=self._interactions_group, interaction=interaction, frame_index=frame_index)
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
            forces[interaction.frame_slice, interaction.particle_indices] += interaction.forces
        return forces

    def __repr__(self) -> str:
        return f"<InteractionsView {len(self)} interaction(s)>"


class HDF5Trajectory(DynamicStructureMethods, TrajectorySource, _HDF5EditableObject):

    def __len__(self) -> int:
        return self.n_frames

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

    _times = HDF5AppendableArray(h5_name="time", title="Time", shape=(), units="picoseconds")
    times = _times.as_numpy()

    _realtimes = HDF5AppendableArray(
        h5_name="realtime", title="Real Time", shape=(), units="seconds"
    )
    realtimes = _realtimes.as_numpy()

    @classmethod
    def load_file(cls, filename: str):
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
    ):
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

        return obj

    @classmethod
    def in_memory(
            cls,
            *,
            title: Optional[str] = None,
            author: Optional[str] = None,
    ):
        with NamedTemporaryFile() as tmp:
            filename = tmp.name + ".h5"
        file = tables.open_file(filename, mode="w", title=title or "", driver="H5FD_CORE",
                                driver_core_backing_store=0)

        obj = cls(file, True)

        obj.author = author
        obj.application = "narupatools"
        obj.convention_version = "1.1"
        obj.narupatools_convention_version = "1.0"
        obj.program = "narupatools"
        obj.conventions = "Pande,NarupaTools"
        obj.program_version = narupatools.__version__


        return obj

    def __init__(self, file: File, writable: bool):
        self._file = file
        self._writable = writable
        self._start_realtime: Optional[float] = None
        self._topology: Optional[HDF5Topology] = None
        self.interactions = InteractionsView(self)

    @property
    def writable(self) -> bool:
        return self._writable

    @property
    def hdf5_group(self) -> Group:
        return self._file.root

    @property
    def n_atoms(self):
        if hasattr(self, "_n_atoms"):
            return self._n_atoms
        self._n_atoms = int(self._positions.shape[1])
        return self._n_atoms

    @property
    def n_frames(self):
        return self.positions.shape[0]

    def save_frame(
            self,
            *,
            positions: np.ndarray,
            time: float,
            potential_energy: float,
            kinetic_energy: float,
            velocities: Optional[np.ndarray] = None,
            forces: Optional[np.ndarray] = None,
    ) -> None:
        """
        Save a trajectory frame.

        :param positions: Atomic positions in nanometers.
        :param velocities: Atomic velocities in nanometers per picoseconds.
        :param time: Time in picoseconds.
        :param potential_energy: Potential energy in kilojoules per mole.
        :param kinetic_energy: Kinetic energy in kilojoules per mole.
        :param forces: Atomic forces in kilojoules per mole per nanometer.
        """
        if not self._writable:
            raise ValueError("Trajectory is read only")

        if not hasattr(self, "_n_atoms"):
            self._n_atoms = len(positions)

        self._positions.append(np.array([positions]))
        if velocities is not None:
            self._velocities.append(np.array([velocities]))
        if forces is not None:
            self._forces.append(np.array([forces]))
        self._times.append(np.array([time]))
        if self._start_realtime is None:
            self._start_realtime = monotonic()
            self._realtimes.append(np.array([0]))
        else:
            self._realtimes.append(np.array([monotonic() - self._start_realtime]))

        self._potential_energies.append(np.array([potential_energy]))
        self._kinetic_energies.append(np.array([kinetic_energy]))

        self._file.flush()

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

    @property
    def topology(self) -> HDF5Topology:
        if self._topology is None:
            if 'topology' in self._file.root:
                self._topology = HDF5Topology.from_string(self._file.root.topology[0].decode("ascii"))
            else:
                raise AttributeError
        return self._topology

    @property
    def masses(self) -> np.ndarray:
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

    def save_to_file(self, filename: str, overwrite_existing=False):
        self._file.copy_file(filename, overwrite_existing)

    @classmethod
    @contextlib.contextmanager
    def record(cls, dynamics: InteractiveSimulationDynamics, *, filename: Optional[str] = None) -> HDF5Trajectory:
        """Record dynamics to a single trajectory, stopping if it is reset."""
        if filename is not None:
            traj = cls.new_file(filename=filename)
        else:
            traj = cls.in_memory()

        has_logged_initial = False
        frame_index = 0

        def log_initial(**kwargs: Any) -> None:
            nonlocal traj
            nonlocal has_logged_initial

            if has_logged_initial:
                return

            log_step(**kwargs)
            traj.save_topology(dynamics.get_frame(fields=everything()))
            has_logged_initial = True

        def log_step(**kwargs: Any) -> None:
            nonlocal traj
            nonlocal frame_index

            traj.save_frame(
                positions=dynamics.positions,
                velocities=dynamics.velocities,
                forces=dynamics.forces,
                kinetic_energy=dynamics.kinetic_energy,
                potential_energy=dynamics.potential_energy,
                time=dynamics.elapsed_time,
            )
            for interaction in dynamics.imd.current_interactions.values():
                traj.interactions.save_interaction(
                    key=interaction.key,
                    interaction=interaction,
                    frame_index=frame_index,
                )

            frame_index += 1

        def log_end_interaction(key: str, **kwargs: Any) -> None:
            nonlocal traj
            nonlocal frame_index

            traj.interactions.end_interaction(key=key, frame_index=frame_index+1)

        dynamics.on_pre_step.add_callback(log_initial, priority=-1000)
        dynamics.on_post_step.add_callback(log_step, priority=-1000)
        dynamics.imd.on_end_interaction.add_callback(log_end_interaction)

        try:
            yield traj
        finally:
            dynamics.on_pre_step.remove_callback(log_initial)
            dynamics.on_post_step.remove_callback(log_step)
            dynamics.imd.on_end_interaction.remove_callback(log_end_interaction)
            traj._writable = False

    def __str__(self):
        s = "<HDF5Trajectory"
        if not self._writable:
            s += " read-only"
        s += f" {len(self.positions)} frame(s)"
        s += f" {self.n_atoms} atom(s)"
        s += ">"
        return s
