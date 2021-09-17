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
from typing import AbstractSet, Any, Dict, List, Union

import numpy as np
from infinite_sets import InfiniteSet, everything
from narupa.trajectory import FrameData
from simtk.openmm import CustomExternalForce, System
from simtk.openmm.app import Simulation

from narupatools.imd import InteractiveSimulationDynamics
from narupatools.imd._feature_setclear import SetAndClearInteractionFeature
from narupatools.physics.typing import (
    ScalarArray,
    Vector3,
    Vector3Array,
    Vector3ArrayLike,
)

from ._converter import (
    get_openmm_masses,
    openmm_context_to_frame,
    openmm_topology_to_frame,
)
from ._serializer import deserialize_simulation
from ..override import override


class OpenMMDynamics(InteractiveSimulationDynamics):
    """Dynamics based on an OpenMM simulation."""

    @override
    def _step_internal(self) -> None:
        with self._simulation_lock:
            self._simulation.step(steps=1)

    @override
    def _reset_internal(self) -> None:
        with self._simulation_lock, BytesIO(self._checkpoint) as bytesio:
            self._simulation.loadCheckpoint(bytesio)
            self._simulation.context.reinitialize(preserveState=True)

    @override
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
    @override
    def imd(self) -> OpenMMIMDFeature:  # noqa: D102
        return self._imd

    @property
    @override
    def timestep(self) -> float:  # noqa: D102
        return self._simulation.integrator.getStepSize()._value

    @property
    @override
    def masses(self) -> ScalarArray:  # noqa: D102
        return self._masses

    @property
    @override
    def positions(self) -> Vector3Array:  # noqa: D102
        with self._simulation_lock:
            state = self._simulation.context.getState(getPositions=True)
        return state.getPositions(asNumpy=True)._value

    @positions.setter
    def positions(self, value: Vector3ArrayLike) -> None:
        with self._simulation_lock:
            self._simulation.context.setPositions(np.asfarray(value))

    @property
    @override
    def velocities(self) -> Vector3Array:  # noqa: D102
        with self._simulation_lock:
            state = self._simulation.context.getState(getVelocities=True)
        return state.getVelocities(asNumpy=True)._value

    @velocities.setter
    def velocities(self, value: Vector3ArrayLike) -> None:
        with self._simulation_lock:
            self._simulation.context.setVelocities(np.asfarray(value))

    @property
    @override
    def forces(self) -> Vector3Array:  # noqa: D102
        with self._simulation_lock:
            state = self._simulation.context.getState(getForces=True)
        return state.getForces(asNumpy=True)._value

    @property
    @override
    def kinetic_energy(self) -> float:  # noqa: D102
        with self._simulation_lock:
            state = self._simulation.context.getState(getEnergy=True)
        return state.getKineticEnergy()._value

    @property
    @override
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


def _get_or_create_imd_force(simulation: Simulation) -> CustomExternalForce:
    """Get a suitable force in an OpenMM simulation for IMD."""
    potential_imd_forces = get_imd_forces_from_system(simulation.system)
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


class OpenMMIMDFeature(SetAndClearInteractionFeature[OpenMMDynamics]):
    """Controls interactions applied to an OpenMM simulation."""

    def __init__(self, dynamics: OpenMMDynamics):
        super().__init__(dynamics)
        self._forces_dirty = True
        self.imd_force = _get_or_create_imd_force(dynamics._simulation)

    @override
    def _on_pre_step(self, **kwargs: Any) -> None:
        super()._on_pre_step(**kwargs)
        self._calculate_and_apply_interactions()
        if self._forces_dirty:
            self.imd_force.updateParametersInContext(self.dynamics._simulation.context)
        self._forces_dirty = False

    @override
    def _set_forces(self, forces: Dict[int, Vector3]) -> None:
        for index, force in forces.items():
            self.imd_force.setParticleParameters(int(index), int(index), force)
            self._forces_dirty = True

    @override
    def _clear_forces(self, indices: InfiniteSet[int] = everything()) -> None:
        all_indices: AbstractSet[int] = (
            set(range(len(self.dynamics.masses))) & indices  # type: ignore
        )
        for index in all_indices:
            self.imd_force.setParticleParameters(int(index), int(index), (0, 0, 0))
            self._forces_dirty = True

    @override
    @property
    def _system_size(self) -> int:
        return len(self.dynamics.masses)


IMD_FORCE_EXPRESSION = "-fx * x - fy * y - fz * z"


def create_imd_force() -> CustomExternalForce:
    """
    Returns an empty OpenMM force to communicate imd forces.

    Each particle in the system has a ``fx``, ``fy``, and ``fz`` parameter to
    provide the arbitrary force components.

    The force needs to be populated to include all the particle in the
    simulation :class:`mm.System`.

    .. seealso: populate_imd_force
    """
    force = CustomExternalForce(IMD_FORCE_EXPRESSION)
    force.addPerParticleParameter("fx")
    force.addPerParticleParameter("fy")
    force.addPerParticleParameter("fz")
    return force


def populate_imd_force(force: CustomExternalForce, system: System) -> None:
    """
    Add all the particles to the iMD force.

    The iMD force must be one generated by :func:`create_imd_force`.

    .. seealso: create_imd_force
    """
    # Attach all the particles to the force object, and set the imd force to 0
    for particle in range(system.getNumParticles()):
        force.addParticle(particle, (0, 0, 0))


def add_imd_force_to_system(system: System) -> CustomExternalForce:
    """
    Generate an OpenMM force that accepts arbitrary forces per particle.

    The force is created, populated, added to the system and returned.

    This is the force that is used to communicate the particle interactions from
    Narupa by :class:`NarupaImdReporter`.

    .. seealso: create_imd_force, populate_imd_force
    """
    force = create_imd_force()
    populate_imd_force(force, system)
    system.addForce(force)
    return force


def get_imd_forces_from_system(system: System) -> List[CustomExternalForce]:
    """
    Find the forces that are compatible with an imd force in a given system.

    A compatible force has the expected energy expression, and contains as
    many particles as the system.

    All the compatible force objects are returned.
    """
    system_num_particles = system.getNumParticles()
    return [
        force
        for force in system.getForces()
        if isinstance(force, CustomExternalForce)
        and force.getEnergyFunction() == IMD_FORCE_EXPRESSION
        and force.getNumParticles() == system_num_particles
    ]
