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

"""Writer for the NarupaTools HDF5 format."""

from __future__ import annotations

import contextlib
import os
from time import monotonic as time_monotonic
from types import TracebackType
from typing import Any, Dict, Optional, Tuple, Type

import numpy as np
import tables
from infinite_sets import everything
from narupa.trajectory import FrameData
from tables import EArray, Float32Atom, Group

import narupatools
from narupatools.imd import Interaction, InteractiveSimulationDynamics

from ._utils import generate_topology


class HDF5Writer:
    """Writer for HDF5 trajectory containing interactive forces."""

    def __init__(
        self,
        filename: str,
        *,
        title: Optional[str] = None,
        author: Optional[str] = None,
        overwrite_existing: bool = False,
    ):
        self._setup(
            filename, title=title, author=author, overwrite_existing=overwrite_existing
        )

    def _setup(
        self,
        filename: str,
        *,
        title: Optional[str] = None,
        author: Optional[str] = None,
        overwrite_existing: bool = False,
    ) -> None:
        if not overwrite_existing and os.path.exists(filename):
            raise FileExistsError(filename)
        self._file = tables.open_file(filename, mode="w", title=title or "")
        self._add_global_attributes()
        if author is not None:
            self._file.root._v_attrs["author"] = str(author)
        self._coordinates: Optional[EArray] = None
        self._forces: Optional[EArray] = None
        self._velocities: Optional[EArray] = None
        self._time: Optional[EArray] = None
        self._realtime: Optional[EArray] = None
        self._kinetic_energy: Optional[EArray] = None
        self._potential_energy: Optional[EArray] = None
        self._interactions: Optional[Group] = None
        self._has_saved_topology: bool = False
        self._start_realtime: Optional[float] = None
        self._per_interaction: Dict[str, Group] = {}
        self._closed = False

    def reopen(
        self,
        filename: str,
        *,
        title: Optional[str] = None,
        author: Optional[str] = None,
    ) -> None:
        """
        Reopen a writer for a different file, discarding the old state.

        This closes the previous file.

        :param filename: Name of the file to write to.
        :param title: Title to add to the trajectory.
        :param author: Name of the author.
        """
        if not self._closed:
            self.close()
        self._setup(filename, title=title, author=author)

    def __enter__(self) -> HDF5Writer:
        return self

    def __exit__(
        self,
        t: Optional[Type[BaseException]],
        value: Optional[BaseException],
        traceback: Optional[TracebackType],
    ) -> None:
        self.close()

    def _add_global_attributes(self) -> None:
        self._file.root._v_attrs["application"] = "narupatools"
        self._file.root._v_attrs["conventionVersion"] = "1.1"
        self._file.root._v_attrs["narupatoolsConventionVersion"] = "1.0"
        self._file.root._v_attrs["program"] = "narupatools"
        self._file.root._v_attrs["conventions"] = "Pande,NarupaTools"
        self._file.root._v_attrs["programVersion"] = narupatools.__version__

    def _create_array(
        self,
        *,
        name: str,
        title: str,
        shape: Tuple[int, ...],
        units: Optional[str] = None,
        parent: Optional[Group] = None,
    ) -> EArray:
        if parent is None:
            parent = self._file.root
        array = self._file.create_earray(
            parent,
            name=name,
            atom=Float32Atom(),
            shape=(0,) + shape,
            title=title,
            filters=tables.Filters(shuffle=True, complib="zlib", complevel=1),
        )
        if units is not None:
            array._v_attrs["units"] = units
        return array

    @property
    def time(self) -> EArray:
        """Array of frame times in picoseconds."""
        if self._time is None:
            self._time = self._create_array(
                name="time", title="Time", shape=(), units="picoseconds"
            )
        return self._time

    @property
    def realtime(self) -> EArray:
        """Array of real times in seconds since the first frame."""
        if self._realtime is None:
            self._realtime = self._create_array(
                name="realtime", title="Real Time", shape=(), units="seconds"
            )
        return self._realtime

    @property
    def coordinates(self) -> EArray:
        """Array of atom coordinates in nanometers."""
        if self._coordinates is None:
            self._coordinates = self._create_array(
                name="coordinates",
                title="Atomic Coordinates",
                shape=(self.n_atoms, 3),
                units="nanometers",
            )
        return self._coordinates

    @property
    def velocities(self) -> EArray:
        """Array of atom velocities in nanometers per picoseconds."""
        if self._velocities is None:
            self._velocities = self._create_array(
                name="velocities",
                title="Atomic Velocities",
                shape=(self.n_atoms, 3),
                units="nanometer/picosecond",
            )
        return self._velocities

    @property
    def forces(self) -> EArray:
        """Array of forces in kilojoules per mole per nanometer."""
        if self._forces is None:
            self._forces = self._create_array(
                name="forces",
                title="Atomic Forces",
                shape=(self.n_atoms, 3),
                units="kilojoules/(mole*nanometer)",
            )
        return self._forces

    @property
    def kinetic_energy(self) -> EArray:
        """Array of kinetic energies in kilojoules per mole."""
        if self._kinetic_energy is None:
            self._kinetic_energy = self._create_array(
                name="kineticEnergy",
                title="Kinetic Energy",
                shape=(),
                units="kilojoules/mole",
            )
        return self._kinetic_energy

    @property
    def potential_energy(self) -> EArray:
        """Array of potential energies in kilojoules per mole."""
        if self._potential_energy is None:
            self._potential_energy = self._create_array(
                name="potentialEnergy",
                title="Potential Energy",
                shape=(),
                units="kilojoules/mole",
            )
        return self._potential_energy

    @property
    def interactions(self) -> Group:
        """Group containing interactions."""
        if self._interactions is None:
            self._interactions = self._file.create_group(
                self._file.root, name="interactions", title="Interactions"
            )
        return self._interactions

    def close(self) -> None:
        """Close the file being written to."""
        for key in list(self._per_interaction.keys()):
            self.end_interaction(
                key=key, frame_index=len(self._coordinates) or 0  # type: ignore
            )
        self._file.close()
        self._closed = True

    def save_topology(self, frame: FrameData) -> None:
        """
        Save the topology of the given frame to the file.

        :param frame: Frame to get topology from.
        """
        raw_string = generate_topology(frame)
        self._file.create_array(
            self._file.root, name="topology", title="Topology", obj=[raw_string]
        )

    def save_frame(
        self,
        *,
        coordinates: np.ndarray,
        time: float,
        potential_energy: float,
        kinetic_energy: float,
        velocities: Optional[np.ndarray] = None,
        forces: Optional[np.ndarray] = None,
    ) -> None:
        """
        Save a trajectory frame.

        :param coordinates: Atomic coordinates in nanometers.
        :param velocities: Atomic velocities in nanometers per picoseconds.
        :param time: Time in picoseconds.
        :param potential_energy: Potential energy in kilojoules per mole.
        :param kinetic_energy: Kinetic energy in kilojoules per mole.
        :param forces: Atomic forces in kilojoules per mole per nanometer.
        """
        if not hasattr(self, "n_atoms"):
            self.n_atoms = len(coordinates)

        self.coordinates.append(np.array([coordinates]))
        if velocities is not None:
            self.velocities.append(np.array([velocities]))
        if forces is not None:
            self.forces.append(np.array([forces]))
        self.time.append(np.array([time]))
        if self._start_realtime is None:
            self._start_realtime = time_monotonic()
            self.realtime.append(np.array([0]))
        else:
            self.realtime.append(np.array([time_monotonic() - self._start_realtime]))

        self.potential_energy.append(np.array([potential_energy]))
        self.kinetic_energy.append(np.array([kinetic_energy]))

        self._file.flush()

    def create_interaction_if_absent(
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
        if key in self._per_interaction:
            return

        hdf_interaction = self._file.create_group(
            self.interactions, name=key, title="Interaction"
        )
        self._per_interaction[key] = hdf_interaction

        hdf_interaction._v_attrs["type"] = interaction.interaction_type
        hdf_interaction._v_attrs["startIndex"] = frame_index

        interaction_size = len(interaction.particle_indices)

        self._file.create_array(
            hdf_interaction,
            name="indices",
            title="Particle Indices",
            obj=list(interaction.particle_indices),
        )

        self._create_array(
            name="position",
            title="Interaction Position",
            shape=(3,),
            units="nanometers",
            parent=hdf_interaction,
        )

        self._create_array(
            name="forces",
            title="Atomic Forces",
            shape=(interaction_size, 3),
            units="kilojoules/(mole*nanometer)",
            parent=hdf_interaction,
        )
        self._create_array(
            name="potentialEnergy",
            title="Potential Energy",
            shape=(),
            units="kilojoules/(mole*nanometer)",
            parent=hdf_interaction,
        )
        self._create_array(
            name="frameIndex",
            title="Frame Index",
            shape=(),
            units=None,
            parent=hdf_interaction,
        )
        self._create_array(
            name="scale", title="Scale", shape=(), units=None, parent=hdf_interaction
        )

        self.save_interaction(
            key=key,
            interaction=interaction,
            frame_index=frame_index,
        )

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
        hdf_interaction = self._per_interaction[key]

        hdf_interaction.potentialEnergy.append([interaction.potential_energy])
        hdf_interaction.forces.append([interaction.forces])
        if hasattr(interaction, "position"):
            hdf_interaction.position.append([interaction.position])  # type: ignore[attr-defined]
        if hasattr(interaction, "scale"):
            hdf_interaction.scale.append([interaction.scale])  # type: ignore[attr-defined]
        hdf_interaction.frameIndex.append([frame_index])

    def end_interaction(self, *, key: str, frame_index: int) -> None:
        """
        Mark an interaction as having finished.

        :param key: Key of the interaction.
        :param frame_index: Current frame of the simulation.
        """
        if key not in self._per_interaction:
            return
        self._per_interaction[key]._v_attrs["endIndex"] = frame_index
        del self._per_interaction[key]


def record_hdf5(
    dynamics: InteractiveSimulationDynamics,
    *,
    filename: str,
    title: Optional[str] = None,
    author: Optional[str] = None,
    write_velocities: bool = True,
    overwrite_existing: bool = False,
) -> HDF5Writer:
    """
    Add an HDF5 Writer to interactive dynamics.

    :param dynamics: Dynamics to attach this writer to.
    :param filename: Filename to write to.
    :param title: Optional title for this trajectory.
    :param author: Optional author to store with the trajectory.
    :param write_velocities: Should this writer write velocities?
    :param overwrite_existing: Should an existing file be overwritten?
    :return: HDF5 writer that is listening to the trajectory.
    """
    fullfilename, extension = os.path.splitext(filename)

    writer = HDF5Writer(
        filename=fullfilename + extension,
        title=title,
        author=author,
        overwrite_existing=overwrite_existing,
    )
    has_logged_initial = False

    def log_initial(**kwargs: Any) -> None:
        nonlocal writer
        nonlocal has_logged_initial

        for interaction in dynamics.imd.current_interactions.values():
            writer.create_interaction_if_absent(
                key=interaction.key,
                interaction=interaction,
                frame_index=dynamics.elapsed_steps,
            )

        if has_logged_initial:
            return

        log_step(**kwargs, write_interactions=False)
        writer.save_topology(dynamics.get_frame(everything()))
        has_logged_initial = True

    def log_step(write_interactions = True, **kwargs: Any) -> None:
        nonlocal writer
        writer.save_frame(
            coordinates=dynamics.positions,
            velocities=dynamics.velocities if write_velocities else None,
            forces=dynamics.forces,
            kinetic_energy=dynamics.kinetic_energy,
            potential_energy=dynamics.potential_energy,
            time=dynamics.elapsed_time,
        )
        if write_interactions:
            for interaction in dynamics.imd.current_interactions.values():
                writer.save_interaction(
                    key=interaction.key,
                    interaction=interaction,
                    frame_index=dynamics.elapsed_steps,
                )

    def log_end_interaction(key: str, **kwargs: Any) -> None:
        nonlocal writer
        writer.end_interaction(key=key, frame_index=dynamics.total_steps)

    def add_callbacks():
        dynamics.on_pre_step.add_callback(log_initial)
        dynamics.on_post_step.add_callback(log_step)
        dynamics.on_reset.add_callback(on_reset)
        dynamics.imd.on_end_interaction.add_callback(log_end_interaction)

    def remove_callbacks():
        dynamics.on_pre_step.remove_callback(log_initial)
        dynamics.on_post_step.remove_callback(log_step)
        dynamics.on_reset.remove_callback(on_reset)
        dynamics.imd.on_end_interaction.remove_callback(log_end_interaction)

    def on_reset(**kwargs: Any) -> None:
        nonlocal writer
        nonlocal has_logged_initial
        writer.close()
        suffix = 2
        while os.path.exists(new_filename := f"{fullfilename}-{suffix}{extension}"):
            suffix += 1
        writer.reopen(filename=new_filename, title=title, author=author)
        has_logged_initial = False
        add_callbacks()

    orig_close = writer.close

    add_callbacks()

    def close_wrapped(**kwargs):
        remove_callbacks()
        orig_close()

    writer.close = close_wrapped

    return writer