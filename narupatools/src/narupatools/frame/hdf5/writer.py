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

from types import TracebackType
from typing import Dict, Optional, Tuple, Type

import numpy as np
import tables
from narupa.imd import ParticleInteraction
from narupa.trajectory import FrameData
from tables import EArray, Float32Atom, Group

import narupatools
from .utils import generate_topology


class HDF5Writer:
    """Writer for HDF5 trajectory containing interactive forces."""

    def __init__(self, filename: str, title: str = ""):
        self._file = tables.open_file(filename, mode="w", title=title)
        self._add_global_attributes()
        self._coordinates: Optional[EArray] = None
        self._forces: Optional[EArray] = None
        self._velocities: Optional[EArray] = None
        self._time: Optional[EArray] = None
        self._kinetic_energy: Optional[EArray] = None
        self._potential_energy: Optional[EArray] = None
        self._interactions: Optional[Group] = None
        self._has_saved_topology: bool = False
        self._per_interaction: Dict[str, Group] = {}

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
        parent: Optional[Group] = None
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
        self._file.close()

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
        coordinates: np.ndarray,
        time: float,
        potentialEnergy: float,
        kineticEnergy: float,
        velocities: Optional[np.ndarray] = None,
        forces: Optional[np.ndarray] = None,
    ) -> None:
        """
        Save a trajectory frame.

        :param coordinates: Atomic coordinates in nanometers.
        :param velocities: Atomic velocities in nanometers per picoseconds.
        :param time: Time in picoseconds.
        :param potentialEnergy: Potential energy in kilojoules per mole.
        :param kineticEnergy: Kinetic energy in kilojoules per mole.
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
        self.potential_energy.append(np.array([potentialEnergy]))
        self.kinetic_energy.append(np.array([kineticEnergy]))

        self._file.flush()

    def create_interaction(
        self,
        *,
        key: str,
        interaction: ParticleInteraction,
        frame_index: int,
        potential_energy: Optional[float] = None,
        forces: Optional[np.ndarray] = None
    ) -> None:
        """
        Save an IMD interaction to the trajectory.

        :param frame_index: Index of the frame in the main trajectory this interaction
                            corresponds to.
        :param key: Unique key of interaction.
        :param interaction: Interaction to save.
        :param potential_energy: Potential energy of interaction in kilojoules per mole.
        :param forces: Forces on affected atoms in kilojoules per mole per nanometer.
        """
        if key in self._per_interaction:
            return

        hdf_interaction = self._file.create_group(
            self.interactions, name=key, title="Interaction"
        )
        self._per_interaction[key] = hdf_interaction

        hdf_interaction._v_attrs["type"] = interaction.interaction_type
        hdf_interaction._v_attrs["startIndex"] = frame_index

        interaction_size = len(interaction.particles)

        self._file.create_array(
            hdf_interaction,
            name="indices",
            title="Particle Indices",
            obj=list(interaction.particles),
        )

        self._create_array(
            name="position",
            title="Interaction Position",
            shape=(3,),
            units="nanometers",
            parent=hdf_interaction,
        )
        if forces is not None:
            self._create_array(
                name="forces",
                title="Atomic Forces",
                shape=(interaction_size, 3),
                units="kilojoules/(mole*nanometer)",
                parent=hdf_interaction,
            )
        if potential_energy is not None:
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
            potential_energy=potential_energy,
            forces=forces,
        )

    def save_interaction(
        self,
        *,
        key: str,
        interaction: ParticleInteraction,
        frame_index: int,
        potential_energy: Optional[float] = None,
        forces: Optional[np.ndarray] = None
    ) -> None:
        """
        Save an IMD interaction to the trajectory.

        :param frame_index: Index of the frame in the main trajectory this interaction
                            corresponds to.
        :param key: Unique key of interaction.
        :param interaction: Interaction to save.
        :param potential_energy: Potential energy of interaction in kilojoules per mole.
        :param forces: Forces on affected atoms in kilojoules per mole per nanometer.
        """
        hdf_interaction = self._per_interaction[key]

        if potential_energy is not None:
            hdf_interaction.potentialEnergy.append([potential_energy])
        if forces is not None:
            hdf_interaction.forces.append([forces])
        hdf_interaction.position.append([interaction.position])
        hdf_interaction.scale.append([interaction.scale])
        hdf_interaction.frameIndex.append([frame_index])

    def end_interaction(self, *, key: str, frame_index: int) -> None:
        """
        Mark an interaction as having finished.

        :param key: Key of the interaction.
        :param frame_index: Current frame of the simulation.
        """
        self._per_interaction[key]._v_attrs["endIndex"] = frame_index
        del self._per_interaction[key]
