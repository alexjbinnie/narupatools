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

"""Simulation dynamics implementation using OpenMM. Currently does not support iMD."""

from __future__ import annotations

import warnings
from io import BytesIO
from os import PathLike
from threading import Lock
from typing import AbstractSet, Any, Dict, Tuple, Union

import numpy as np
from infinite_sets import InfiniteSet, everything
from narupa.imd import ParticleInteraction
from narupa.openmm.imd import add_imd_force_to_system, get_imd_forces_from_system
from narupa.trajectory import FrameData
from simtk.openmm import CustomExternalForce
from simtk.openmm.app import Simulation

from narupatools.core.broadcastable import Broadcastable
from narupatools.imd import (
    Interaction,
    InteractiveSimulationDynamics,
    calculate_imd_force,
)
from narupatools.imd.feature_setclear import SetAndClearInteractionFeature
from narupatools.physics.typing import (
    ScalarArray,
    Vector3,
    Vector3Array,
    Vector3ArrayLike,
)
from .converter import (
    get_openmm_masses,
    openmm_context_to_frame,
    openmm_topology_to_frame,
)
from .serializer import deserialize_simulation


class OpenMMDynamics(InteractiveSimulationDynamics, Broadcastable):
    """Dynamics based on an OpenMM simulation."""

    def _step_internal(self) -> None:
        with self._simulation_lock:
            self._simulation.step(steps=1)

    def _reset_internal(self) -> None:
        with self._simulation_lock, BytesIO(self._checkpoint) as bytesio:
            self._simulation.loadCheckpoint(bytesio)
            self._simulation.context.reinitialize(preserveState=True)

    def _get_frame(self, fields: InfiniteSet[str]) -> FrameData:
        frame = FrameData()
        with self._simulation_lock:
            openmm_context_to_frame(
                self._simulation.context, fields=fields, existing=frame
            )
            openmm_topology_to_frame(
                self._simulation.topology, fields=fields, existing=frame
            )
        return frame

    def __init__(self, simulation: Simulation, playback_interval: float = 0.0):
        super().__init__(playback_interval=playback_interval)
        self._simulation = simulation
        self._masses = get_openmm_masses(simulation)

        self._imd = OpenMMIMDFeature(self)
        self._simulation_lock = Lock()

        with BytesIO() as bytesio:
            self._simulation.saveCheckpoint(bytesio)
            self._checkpoint = bytesio.getvalue()

    @property
    def imd(self) -> OpenMMIMDFeature:  # noqa: D102
        return self._imd

    @property
    def timestep(self) -> float:  # noqa: D102
        return self._simulation.integrator.getStepSize()._value

    @property
    def masses(self) -> ScalarArray:  # noqa: D102
        return self._masses

    @property
    def positions(self) -> Vector3Array:  # noqa: D102
        with self._simulation_lock:
            state = self._simulation.context.getState(getPositions=True)
        return state.getPositions(asNumpy=True)._value

    @positions.setter
    def positions(self, value: Vector3ArrayLike) -> None:
        with self._simulation_lock:
            self._simulation.context.setPositions(np.asfarray(value))

    @property
    def velocities(self) -> Vector3Array:  # noqa: D102
        with self._simulation_lock:
            state = self._simulation.context.getState(getVelocities=True)
        return state.getVelocities(asNumpy=True)._value

    @velocities.setter
    def velocities(self, value: Vector3ArrayLike) -> None:
        with self._simulation_lock:
            self._simulation.context.setVelocities(np.asfarray(value))

    @property
    def forces(self) -> Vector3Array:  # noqa: D102
        with self._simulation_lock:
            state = self._simulation.context.getState(getForces=True)
        return state.getForces(asNumpy=True)._value

    @property
    def kinetic_energy(self) -> float:  # noqa: D102
        with self._simulation_lock:
            state = self._simulation.context.getState(getEnergy=True)
        return state.getKineticEnergy()._value

    @property
    def potential_energy(self) -> float:  # noqa: D102
        with self._simulation_lock:
            state = self._simulation.context.getState(getEnergy=True)
        return state.getPotentialEnergy()._value

    @staticmethod
    def from_xml_file(path: Union[str, bytes, PathLike], /) -> OpenMMDynamics:
        """
        Create OpenMM dynamics from an XML file path.

        :param path: Path to an XML file of a serialized OpenMM simulation.
        :return: A dynamics object that can be broadcast on a Narupa session.
        """
        with open(path) as infile:
            return OpenMMDynamics.from_xml_string(infile.read())

    @staticmethod
    def from_xml_string(string: str, /) -> OpenMMDynamics:
        """
        Create OpenMM dynamics from the contents of an XML file.

        :param string: Contents of an XML file of a serialized OpenMM simulation.
        :return: A dynamics object that can be broadcast on a Narupa session.
        """
        simulation = deserialize_simulation(string)
        return OpenMMDynamics.from_simulation(simulation)

    @staticmethod
    def from_simulation(simulation: Simulation, /) -> OpenMMDynamics:
        """
        Create OpenMM dynamics from a simulation.

        :param simulation: OpenMM simulation to wrap.
        :return: Dynamics object that can be broadcast on a Narupa session.
        """
        return OpenMMDynamics(simulation)


class OpenMMInteraction(Interaction[OpenMMDynamics]):
    """IMD Interaction applied to an OpenMM simulation."""

    def update_energy_and_forces(self) -> Tuple[Vector3Array, float]:  # noqa: D102
        return calculate_imd_force(
            interaction=self.interaction,
            positions=self._dynamics.positions[self.particle_indices],
            masses=self._dynamics.masses[self.particle_indices],
        )

    def get_positions(self) -> Vector3Array:  # noqa: D102
        return self._dynamics.positions[self.particle_indices]  # type: ignore


def _get_or_create_imd_force(simulation: Simulation) -> CustomExternalForce:
    """Get a suitable force in an OpenMM simulation for IMD."""
    potential_imd_forces = get_imd_forces_from_system(simulation.system)  # type: ignore
    if len(potential_imd_forces) == 0:
        force = add_imd_force_to_system(simulation.system)
        simulation.context.reinitialize(True)
        return force
    if len(potential_imd_forces) > 1:
        warnings.warn(
            "Simulation contains multiple forces suitable for IMD. "
            "The last one will be used."
        )
    return potential_imd_forces[-1]


class OpenMMIMDFeature(
    SetAndClearInteractionFeature[OpenMMDynamics, OpenMMInteraction]
):
    """Controls interactions applied to an OpenMM simulation."""

    def __init__(self, dynamics: OpenMMDynamics):
        super().__init__(dynamics)
        self._forces_dirty = True
        self.imd_force = _get_or_create_imd_force(dynamics._simulation)

    def _on_pre_step(self, **kwargs: Any) -> None:
        super()._on_pre_step(**kwargs)
        self._calculate_and_apply_interactions()
        if self._forces_dirty:
            self.imd_force.updateParametersInContext(self.dynamics._simulation.context)
        self._forces_dirty = False

    def _set_forces(self, forces: Dict[int, Vector3]) -> None:
        for index, force in forces.items():
            self.imd_force.setParticleParameters(int(index), int(index), force)
            self._forces_dirty = True

    def _clear_forces(self, indices: InfiniteSet[int] = everything()) -> None:
        all_indices: AbstractSet[int] = (
            set(range(len(self.dynamics.masses))) & indices  # type: ignore
        )
        for index in all_indices:
            self.imd_force.setParticleParameters(int(index), int(index), (0, 0, 0))
            self._forces_dirty = True

    def _make_interaction(
        self, key: str, interaction: ParticleInteraction, start_time: float
    ) -> OpenMMInteraction:
        return OpenMMInteraction(
            dynamics=self._dynamics,
            key=key,
            interaction=interaction,
            start_time=start_time,
        )

    @property
    def _system_size(self) -> int:
        return len(self.dynamics.masses)
