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
from contextlib import contextmanager
from os import PathLike
from threading import Lock
from typing import AbstractSet, Any, Dict, Generator, List, Optional, Union

import numpy as np
from infinite_sets import InfiniteSet, everything
from narupa.trajectory import FrameData
from openmm import CustomExternalForce, System
from openmm.app import Simulation

from narupatools.frame import (
    DynamicStructureMethods,
    ParticlePositions,
    ParticleVelocities,
)
from narupatools.imd import InteractiveSimulationDynamics, SetAndClearInteractionFeature
from narupatools.override import override
from narupatools.physics.typing import (
    ScalarArray,
    Vector3,
    Vector3Array,
    Vector3ArrayLike,
)

from ._converter import get_openmm_masses
from ._serializer import deserialize_simulation
from ._simulation import OpenMMSimulation


class OpenMMDynamics(InteractiveSimulationDynamics, DynamicStructureMethods):
    """Dynamics based on an OpenMM simulation."""

    @override(InteractiveSimulationDynamics._step_internal)
    def _step_internal(self) -> None:
        with self._simulation_lock:
            self._simulation.run(steps=1)

    @override(InteractiveSimulationDynamics._reset_internal)
    def _reset_internal(self) -> None:
        with self._simulation_lock:
            self._simulation.load_checkpoint(self._checkpoint)
            self._simulation.context.reinitialize(preserveState=True)

    @override(InteractiveSimulationDynamics._get_frame)
    def _get_frame(
        self, fields: InfiniteSet[str], existing: Optional[FrameData] = None
    ) -> FrameData:
        return self._simulation.get_frame(fields=fields, existing=existing)

    def __init__(
        self,
        simulation: Union[Simulation, OpenMMSimulation],
        playback_interval: float = 0.0,
    ):
        super().__init__(playback_interval=playback_interval)
        if isinstance(simulation, Simulation):
            simulation = OpenMMSimulation.from_simulation(simulation)
        self._simulation = simulation
        self._masses = get_openmm_masses(self._simulation.system)

        self._imd = OpenMMIMDFeature(self)
        self._simulation_lock = Lock()

        self._checkpoint = self._simulation.create_checkpoint()

    def set_reset_state(self) -> None:
        """Set the current state as the state that will be returned to on reset."""
        self._checkpoint = self._simulation.create_checkpoint()

    @contextmanager
    def modify_simulation(self) -> Generator[OpenMMSimulation, None, None]:
        """Context manager which allows modification of the underlying simulation."""
        with self._simulation.modify():
            yield self._simulation
        self._checkpoint = self._simulation.create_checkpoint()

    @property
    def simulation(self) -> OpenMMSimulation:
        """
        Underlying OpenMM Simulation.

        Generally, this should not be modified directly, as it may lead to unexpected issues.
        """
        return self._simulation

    @property
    @override(InteractiveSimulationDynamics.imd)
    def imd(self) -> OpenMMIMDFeature:  # noqa: D102
        return self._imd

    @property
    @override(InteractiveSimulationDynamics.timestep)
    def timestep(self) -> float:  # noqa: D102
        return self._simulation.integrator.getStepSize()._value

    @property
    @override(InteractiveSimulationDynamics.masses)
    def masses(self) -> ScalarArray:  # noqa: D102
        return self._masses

    @masses.setter
    def masses(self, value: ScalarArray) -> None:
        with self._simulation_lock:
            system = self._simulation.system
            for i, m in enumerate(value):
                system.setParticleMass(i, m)
            self._simulation.context.reinitialize(preserveState=True)

    @property
    @override(InteractiveSimulationDynamics.positions)
    def positions(self) -> Vector3Array:  # noqa: D102
        with self._simulation_lock:
            state = self._simulation.context.getState(getPositions=True)
        return state.getPositions(asNumpy=True)._value

    @positions.setter
    def positions(self, value: Vector3ArrayLike) -> None:
        with self._simulation_lock:
            self._simulation.context.setPositions(np.asfarray(value))
        for interaction in self._imd.current_interactions.values():
            interaction.mark_positions_dirty()
        self._on_fields_changed.invoke(fields={ParticlePositions})

    @property  # type: ignore
    @override(InteractiveSimulationDynamics.velocities)
    def velocities(self) -> Vector3Array:  # type: ignore  # noqa: D102
        with self._simulation_lock:
            state = self._simulation.context.getState(getVelocities=True)
        return state.getVelocities(asNumpy=True)._value

    @velocities.setter
    def velocities(self, value: Vector3ArrayLike) -> None:
        with self._simulation_lock:
            self._simulation.context.setVelocities(np.asfarray(value))
        for interaction in self._imd.current_interactions.values():
            interaction.mark_velocities_dirty()
        self._on_fields_changed.invoke(fields={ParticleVelocities})

    @property
    @override(InteractiveSimulationDynamics.forces)
    def forces(self) -> Vector3Array:  # noqa: D102
        with self._simulation_lock:
            state = self._simulation.context.getState(getForces=True)
        return state.getForces(asNumpy=True)._value

    @property
    @override(InteractiveSimulationDynamics.kinetic_energy)
    def kinetic_energy(self) -> float:  # noqa: D102
        with self._simulation_lock:
            state = self._simulation.context.getState(getEnergy=True)
        return state.getKineticEnergy()._value

    @property
    @override(InteractiveSimulationDynamics.potential_energy)
    def potential_energy(self) -> float:  # noqa: D102
        with self._simulation_lock:
            state = self._simulation.context.getState(getEnergy=True)
        return state.getPotentialEnergy()._value

    @staticmethod
    def from_xml_file(
        path: Union[str, bytes, PathLike], /, *, platform: Optional[str] = None
    ) -> OpenMMDynamics:
        """
        Create OpenMM dynamics from an XML file path.

        :param path: Path to an XML file of a serialized OpenMM simulation.
        :return: A dynamics object that can be broadcast on a Narupa session.
        """
        with open(path) as infile:
            return OpenMMDynamics.from_xml_string(infile.read(), platform=platform)

    def minimize(self, tolerance: float, max_iterations: Optional[int] = None) -> None:
        """Minimize the system to a given tolerance."""
        with self._simulation_lock:
            if not max_iterations:
                max_iterations = 0
            self._simulation.minimize(tolerance, max_iterations)
            self._on_fields_changed.invoke(fields={ParticlePositions})

    @staticmethod
    def from_xml_string(
        string: str, /, *, platform: Optional[str] = None
    ) -> OpenMMDynamics:
        """
        Create OpenMM dynamics from the contents of an XML file.

        :param string: Contents of an XML file of a serialized OpenMM simulation.
        :return: A dynamics object that can be broadcast on a Narupa session.
        """
        simulation = deserialize_simulation(string, platform=platform)
        return OpenMMDynamics.from_simulation(simulation)

    @staticmethod
    def from_simulation(
        simulation: Union[OpenMMSimulation, Simulation], /
    ) -> OpenMMDynamics:
        """
        Create OpenMM dynamics from a simulation.

        :param simulation: OpenMM simulation to wrap.
        :return: Dynamics object that can be broadcast on a Narupa session.
        """
        if isinstance(simulation, Simulation):
            simulation = OpenMMSimulation.from_simulation(simulation)
        return OpenMMDynamics(simulation)

    @classmethod
    def _create_from_object(cls, obj: Any) -> OpenMMDynamics:
        if isinstance(obj, Simulation):
            return OpenMMDynamics.from_simulation(obj)
        raise NotImplementedError


def _get_or_create_imd_force(simulation: OpenMMSimulation) -> CustomExternalForce:
    """Get a suitable force in an OpenMM simulation for IMD."""
    potential_imd_forces = get_imd_forces_from_system(simulation.system)
    if len(potential_imd_forces) == 0:
        with simulation.modify():
            force = add_imd_force_to_system(simulation.system)
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

    @override(SetAndClearInteractionFeature._on_pre_step)
    def _on_pre_step(self, **kwargs: Any) -> None:
        super()._on_pre_step(**kwargs)
        self._calculate_and_apply_interactions()
        if self._forces_dirty:
            self.imd_force.updateParametersInContext(self.dynamics._simulation.context)
        self._forces_dirty = False

    @override(SetAndClearInteractionFeature._set_forces)
    def _set_forces(self, forces: Dict[int, Vector3]) -> None:
        for index, force in forces.items():
            self.imd_force.setParticleParameters(int(index), int(index), force)  # type: ignore
            self._forces_dirty = True

    @override(SetAndClearInteractionFeature._clear_forces)
    def _clear_forces(self, indices: InfiniteSet[int] = everything()) -> None:
        all_indices: AbstractSet[int] = (
            set(range(len(self.dynamics.masses))) & indices  # type: ignore
        )
        for index in all_indices:
            self.imd_force.setParticleParameters(int(index), int(index), (0, 0, 0))
            self._forces_dirty = True

    @override(SetAndClearInteractionFeature._system_size)
    @property
    def _system_size(self) -> int:
        return len(self.dynamics.masses)

    def _on_post_step(self, **kwargs: Any) -> None:
        for interaction in self.current_interactions.values():
            interaction.mark_positions_dirty()
            interaction.mark_velocities_dirty()
        super()._on_post_step(**kwargs)


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
