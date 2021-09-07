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

"""Simulation dynamics implementation using ASE."""

from __future__ import annotations

from threading import Lock
from typing import Dict, Generic, TypeVar

import numpy as np
import numpy.typing as npt
from ase.atoms import Atoms
from ase.md import Langevin, VelocityVerlet
from ase.md.md import MolecularDynamics
from infinite_sets import InfiniteSet
from narupa.trajectory import FrameData

from narupatools.core.units import UnitsNarupa
from narupatools.imd import Interaction, InteractiveSimulationDynamics
from narupatools.imd._feature import InteractionFeature
from narupatools.imd.interactions._interactiondata import InteractionData
from narupatools.physics.quaternion import quaternion
from narupatools.physics.typing import ScalarArray, Vector3Array, Vector3ArrayLike
from ._rotations import get_angular_velocities

from ..core.dynamics import SimulationRotationProperties
from ._converter import ase_atoms_to_frame
from ._rotational_velocity_verlet import (
    get_angular_momenta,
    get_principal_moments,
    get_rotations,
    set_angular_momenta,
    set_principal_moments,
    set_rotations,
)
from ._system import ASESystem
from ._units import UnitsASE
from .calculators import NullCalculator
from .constraints import InteractionConstraint

TIntegrator = TypeVar("TIntegrator", bound=MolecularDynamics)

_NarupaToASE = UnitsNarupa >> UnitsASE
_ASEToNarupa = UnitsASE >> UnitsNarupa


class ASEDynamics(
    InteractiveSimulationDynamics, SimulationRotationProperties, Generic[TIntegrator]
):
    """
    Run dynamics using an ASE `MolecularDynamics` object.

    This allows the simulation to be run at a set playback rate, as well as exposing
    standard properties such as time step and elapsed steps.
    """

    def __init__(
        self,
        dynamics: TIntegrator,
        *,
        playback_interval: float = 0.0,
    ):
        """
        Create a `SimulationDynamics` wrapper around an ASE `MolecularDynamics` object.

        :param dynamics: ASE `MolecularDynamics` object to be simulated.
        :param playback_interval: Time to wait between steps in seconds. Defaults to 0.0
                                  (run as fast as possible).
        """
        super().__init__(playback_interval=playback_interval)
        self._dynamics = dynamics
        self._initial_positions = self.atoms.get_positions()
        self._initial_momenta = self.atoms.get_momenta()
        self._initial_box = self.atoms.get_cell()
        self._imd = ASEIMDFeature(self)
        self._atom_lock = Lock()

    @staticmethod
    def from_ase_dynamics(dynamics: TIntegrator) -> ASEDynamics[TIntegrator]:
        """
        Create ASE dynamics from an existing dynamics object.

        This object implements the SimulationDynamics API, and hence can be played at a
        specific playback rate. It also provides callbacks for when the simulation
        advances a step or is reset.

        :param dynamics: Any ASE dynamics derived from MolecularDynamics.
        :return: An ASEDynamics object which wraps the provided dynamics.
        """
        return ASEDynamics(dynamics)

    @staticmethod
    def create_langevin(
        atoms: Atoms,
        *,
        friction: float = 1e-2,
        temperature: float = 300,
        timestep: float = 1,
    ) -> ASEDynamics[Langevin]:
        """
        Create ASE Langevin dynamics with the provided ASE atoms object.

        :param atoms: An ASE atoms object to simulate.
        :param friction: Friction of the Langevin integrator in inverse picoseconds.
        :param temperature: Temperature of the Langevin integrator in kelvin.
        :param timestep: Timestep of the Langevin integrator in picoseconds.
        :return: An ASEDynamics object which wraps a Langevin integrator running on the
                 specified system.
        """
        if atoms.get_calculator() is None:
            atoms.set_calculator(NullCalculator())
        dynamics = Langevin(
            atoms,
            timestep=timestep * _NarupaToASE.time,
            temperature_K=temperature,
            friction=friction / _NarupaToASE.time,
            fixcm=False,
        )
        return ASEDynamics.from_ase_dynamics(dynamics)

    @staticmethod
    def create_velocity_verlet(
        atoms: Atoms, timestep: float = 1
    ) -> ASEDynamics[VelocityVerlet]:
        """
        Create ASE velocity verlet dynamics for the provided ASE atoms object.

        :param atoms: An ASE atoms object to simulate.
        :param timestep: Timestep of the Velocity Verlet integrator in picoseconds.
        :return: An ASEDynamics object which wraps a Velocity Verlet integrator running
                 on the specified system.
        """
        if atoms.get_calculator() is None:
            atoms.set_calculator(NullCalculator())
        dynamics = VelocityVerlet(atoms, timestep=timestep * _NarupaToASE.time)
        return ASEDynamics.from_ase_dynamics(dynamics)

    @property
    def imd(self) -> ASEIMDFeature:  # noqa:D102
        return self._imd

    @property
    def atoms(self) -> Atoms:
        """
        Get the wrapped ASE `Atoms` object.

        Direct modification of this may have side effects.
        """
        return self.molecular_dynamics.atoms

    @property
    def molecular_dynamics(self) -> TIntegrator:
        """
        Get the wrapped ASE `MolecularDynamics`` object.

        Direct modification of this may have side effects.
        """
        return self._dynamics

    def _step_internal(self) -> None:
        with self._atom_lock:
            self._dynamics.run(1)

    def _reset_internal(self) -> None:
        with self._atom_lock:
            self.atoms.set_positions(self._initial_positions)
            self.atoms.set_momenta(self._initial_momenta)
            self.atoms.set_cell(self._initial_box)

    @property
    def timestep(self) -> float:  # noqa: D102
        return self.molecular_dynamics.dt * _ASEToNarupa.time

    @timestep.setter
    def timestep(self, value: float) -> None:
        self.molecular_dynamics.dt = value * _NarupaToASE.time

    def _get_frame(self, fields: InfiniteSet[str]) -> FrameData:
        frame = FrameData()
        with self._atom_lock:
            ase_atoms_to_frame(self.atoms, fields=fields, frame=frame)
        return frame

    @property
    def positions(self) -> Vector3Array:  # noqa: D102
        return self.atoms.positions * _ASEToNarupa.length  # type: ignore

    @positions.setter
    def positions(self, value: Vector3ArrayLike) -> None:
        self.atoms.set_positions(np.asfarray(value) * _NarupaToASE.length)

    @property
    def velocities(self) -> Vector3Array:  # noqa: D102
        return self.atoms.get_velocities() * _ASEToNarupa.velocity  # type: ignore

    @velocities.setter
    def velocities(self, value: Vector3ArrayLike) -> None:
        self.atoms.set_velocities(np.asfarray(value) * _NarupaToASE.velocity)

    @property
    def forces(self) -> Vector3Array:  # noqa: D102
        return self.atoms.get_forces() * _ASEToNarupa.force  # type: ignore

    @property
    def masses(self) -> ScalarArray:  # noqa: D102
        return self.atoms.get_masses() * _ASEToNarupa.mass  # type: ignore

    @property
    def kinetic_energy(self) -> float:  # noqa: D102
        return self.atoms.get_kinetic_energy() * _ASEToNarupa.energy

    @property
    def potential_energy(self) -> float:  # noqa: D102
        return self.atoms.get_potential_energy() * _ASEToNarupa.energy

    @property
    def orientations(self) -> npt.NDArray[quaternion]:  # noqa: D102
        return get_rotations(self.atoms)

    @orientations.setter
    def orientations(self, value: npt.NDArray[quaternion]) -> None:
        set_rotations(self.atoms, value)

    @property
    def angular_momenta(self) -> Vector3Array:  # noqa: D102
        return get_angular_momenta(self.atoms) * _ASEToNarupa.angular_momentum

    @angular_momenta.setter
    def angular_momenta(self, value: Vector3ArrayLike) -> None:
        set_angular_momenta(
            self.atoms, np.asfarray(value) * _NarupaToASE.angular_momentum
        )

    @property
    def angular_velocities(self) -> Vector3Array:  # noqa: D102
        """Angular velocity of each particle abouts its center of mass, in radians per picoseconds."""
        return get_angular_velocities(self.atoms) * _ASEToNarupa.angular_velocity

    @property
    def moments_of_inertia(self) -> Vector3Array:  # noqa: D102
        return get_principal_moments(self.atoms) * _ASEToNarupa.moment_inertia

    @moments_of_inertia.setter
    def moments_of_inertia(self, value: Vector3ArrayLike) -> None:
        set_principal_moments(
            self.atoms, np.asfarray(value) * _NarupaToASE.moment_inertia
        )

    @property
    def torques(self) -> Vector3Array: ## noqa: D102
        return self.atoms.get_torques() * _NarupaToASE.torque


class ASEIMDFeature(InteractionFeature[ASEDynamics]):
    """Interactive Molecular Dynamics manager for ASE dynamics."""

    def __init__(self, dynamics: ASEDynamics):
        super().__init__(dynamics)
        self.constraints: Dict[str, InteractionConstraint] = {}

    @property
    def _system_size(self) -> int:
        return len(self.dynamics.atoms)

    def create_interaction(  # noqa: D102
        self, *, key: str, interaction: InteractionData, start_time: float
    ) -> Interaction:
        # The interaction does not store a reference to the dynamics itself.
        instance = Interaction.create(
            key=key,
            interaction=interaction,
            start_time=start_time,
            dynamics=ASESystem(self.dynamics.atoms),
        )
        constraint = InteractionConstraint(interaction=instance)
        self.constraints[key] = constraint
        self.dynamics.atoms.constraints.append(constraint)
        return instance

    def remove_interaction(self, key: str) -> Interaction:  # noqa: D102
        instance = super().remove_interaction(key)
        constraint = self.constraints[key]
        self.dynamics.atoms.constraints.remove(constraint)
        return instance
