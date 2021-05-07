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
from typing import Generic, Optional, TypeVar

import numpy as np
from MDAnalysis import Universe
from ase.atoms import Atoms
from ase.md import Langevin, VelocityVerlet
from ase.md.md import MolecularDynamics
from infinite_sets import InfiniteSet
from narupa.trajectory import FrameData

from narupatools.ase import NullCalculator
from narupatools.core.units import UnitsNarupa
from narupatools.imd import InteractiveSimulationDynamics
from narupatools.mdanalysis import mdanalysis_universe_to_frame
from narupatools.physics.typing import ScalarArray, Vector3Array, Vector3ArrayLike
from .converter import ase_atoms_to_frame
from .imd import ASEIMDFeature
from .units import UnitsASE

TIntegrator = TypeVar("TIntegrator", bound=MolecularDynamics)

_NarupaToASE = UnitsNarupa >> UnitsASE
_ASEToNarupa = UnitsASE >> UnitsNarupa


class ASEDynamics(InteractiveSimulationDynamics, Generic[TIntegrator]):
    """
    Run dynamics using an ASE `MolecularDynamics` object.

    This allows the simulation to be run at a set playback rate, as well as exposing
    standard properties such as time step and elapsed steps.
    """

    def __init__(
        self,
        dynamics: TIntegrator,
        *,
        universe: Optional[Universe] = None,
        playback_interval: float = 0.0,
    ):
        """
        Create a `SimulationDynamics` wrapper around an ASE `MolecularDynamics` object.

        :param dynamics: ASE `MolecularDynamics` object to be simulated.
        :param universe: Optional MDAnalysis universe to provide additional topology.
        :param playback_interval: Time to wait between steps in seconds. Defaults to 0.0
                                  (run as fast as possible).
        """
        super().__init__(playback_interval=playback_interval)
        self._dynamics = dynamics
        self._universe = universe
        self._initial_positions = self.atoms.get_positions()
        self._initial_momenta = self.atoms.get_momenta()
        self._initial_box = self.atoms.get_cell()
        self._imd = ASEIMDFeature(self)
        self._atom_lock = Lock()

    @staticmethod
    def from_ase_dynamics(
        dynamics: TIntegrator, *, universe: Optional[Universe] = None
    ) -> ASEDynamics[TIntegrator]:
        """
        Create ASE dynamics from an existing dynamics object.

        This object implements the SimulationDynamics API, and hence can be played at a
        specific playback rate. It also provides callbacks for when the simulation
        advances a step or is reset.

        :param dynamics: Any ASE dynamics derived from MolecularDynamics.
        :param universe: Optional MDAnalysis universe with additional system
                         information.
        :return: An ASEDynamics object which wraps the provided dynamics.
        """
        return ASEDynamics(dynamics, universe=universe)

    @staticmethod
    def create_langevin(
        atoms: Atoms,
        *,
        universe: Optional[Universe] = None,
        friction: float = 1e-2,
        temperature: float = 300,
        timestep: float = 1,
    ) -> ASEDynamics[Langevin]:
        """
        Create ASE Langevin dynamics with the provided ASE atoms object.

        :param atoms: An ASE atoms object to simulate.
        :param universe: Optional MDAnalysis universe with additional system
                         information.
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
        return ASEDynamics.from_ase_dynamics(dynamics, universe=universe)

    @staticmethod
    def create_velocity_verlet(
        atoms: Atoms, *, universe: Optional[Universe] = None, timestep: float = 1
    ) -> ASEDynamics[VelocityVerlet]:
        """
        Create ASE velocity verlet dynamics for the provided ASE atoms object.

        :param atoms: An ASE atoms object to simulate.
        :param universe: Optional MDAnalysis universe with additional system
                         information.
        :param timestep: Timestep of the Velocity Verlet integrator in picoseconds.
        :return: An ASEDynamics object which wraps a Velocity Verlet integrator running
                 on the specified system.
        """
        if atoms.get_calculator() is None:
            atoms.set_calculator(NullCalculator())
        dynamics = VelocityVerlet(atoms, timestep=timestep * _NarupaToASE.time)
        return ASEDynamics.from_ase_dynamics(dynamics, universe=universe)

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
        if self._universe:
            with self._atom_lock:
                ase_atoms_to_frame(self.atoms, fields=fields, frame=frame)
            added_fields = set(frame.arrays.keys()) | set(frame.values.keys())
            mdanalysis_universe_to_frame(
                self._universe, fields=fields - added_fields, frame=frame
            )
            return frame
        else:
            with self._atom_lock:
                ase_atoms_to_frame(self.atoms, fields=fields, frame=frame)
        return frame

    @property
    def positions(self) -> Vector3Array:  # noqa: D102
        return self.atoms.positions * _ASEToNarupa.length

    @positions.setter
    def positions(self, value: Vector3ArrayLike) -> None:
        self.atoms.set_positions(np.asfarray(value) * _NarupaToASE.length)

    @property
    def velocities(self) -> Vector3Array:  # noqa: D102
        return self.atoms.get_velocities() * _ASEToNarupa.velocity

    @velocities.setter
    def velocities(self, value: Vector3ArrayLike) -> None:
        self.atoms.set_velocities(np.asfarray(value) * _NarupaToASE.velocity)

    @property
    def forces(self) -> Vector3Array:  # noqa: D102
        return self.atoms.get_forces() * _ASEToNarupa.force

    @property
    def masses(self) -> ScalarArray:  # noqa: D102
        return self.atoms.get_masses() * _ASEToNarupa.mass

    @property
    def kinetic_energy(self) -> float:  # noqa: D102
        return self.atoms.get_kinetic_energy() * _ASEToNarupa.energy

    @property
    def potential_energy(self) -> float:  # noqa: D102
        return self.atoms.get_potential_energy() * _ASEToNarupa.energy
