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

"""Molecular dynamics run directly in LAMMPS."""

from __future__ import annotations

from typing import Optional, Tuple, List

import numpy as np
from infinite_sets import InfiniteSet
from narupa.trajectory import FrameData

from narupatools.imd import InteractionFeature, InteractiveSimulationDynamics
from narupatools.physics.typing import ScalarArray, Vector3Array, Vector3ArrayLike
from narupatools.physics.units import UnitsNarupa, UnitSystem

from ._simulation import LAMMPSSimulation


class LAMMPSDynamics(InteractiveSimulationDynamics):
    """Molecular dynamics implementation using LAMMPS directly."""

    def __init__(
        self,
        simulation: LAMMPSSimulation,
        *,
        playback_interval: float = 0.0,
        units: Optional[UnitSystem] = None,
    ):
        """
        Create a new LAMMPS dynamics for a given simulation.

        :param simulation: LAMMPS simulation to run.
        :param playback_interval: Interval at which to run dynamics, in seconds.
        """
        super().__init__(playback_interval=playback_interval)
        self._simulation = simulation
        self._imd = LAMMPSIMDFeature(self)

        self._lammps_to_narupa = self._simulation.unit_system >> UnitsNarupa
        self._narupa_to_lammps = UnitsNarupa >> self._simulation.unit_system

        simulation.setup_imd(self._imd.calculate_imd)

    @property
    def imd(self) -> LAMMPSIMDFeature:  # noqa:D102
        return self._imd

    @property
    def simulation(self) -> LAMMPSSimulation:
        """Underlying LAMMPS simulation."""
        return self._simulation

    def _step_internal(self) -> None:
        self._simulation.run(1)

    def _reset_internal(self) -> None:
        # Allow subclasses to override behaviour on reset
        pass

    def _get_frame(
        self, fields: InfiniteSet[str], existing: Optional[FrameData] = None
    ) -> FrameData:
        return self._simulation.get_frame(fields=fields, existing=existing)

    @property
    def timestep(self) -> float:  # noqa: D102
        return self._simulation.timestep * self._lammps_to_narupa.time

    @property
    def positions(self) -> Vector3Array:  # noqa: D102
        return self._simulation.positions * self._lammps_to_narupa.length

    @positions.setter
    def positions(self, value: Vector3ArrayLike) -> None:
        self._simulation.positions = np.asfarray(value) * self._narupa_to_lammps.length

    @property
    def velocities(self) -> Vector3Array:  # noqa: D102
        return self._simulation.velocities * self._lammps_to_narupa.velocity

    @velocities.setter
    def velocities(self, value: Vector3ArrayLike) -> None:
        self._simulation.velocities = (
            np.asfarray(value) * self._narupa_to_lammps.velocity
        )

    @property
    def forces(self) -> Vector3Array:  # noqa: D102
        return self._simulation.forces * self._lammps_to_narupa.force

    @property
    def masses(self) -> ScalarArray:  # noqa: D102
        return self._simulation.masses * self._lammps_to_narupa.mass

    @property
    def kinetic_energy(self) -> float:  # noqa: D102
        return self._simulation.kinetic_energy * self._lammps_to_narupa.energy

    @property
    def potential_energy(self) -> float:  # noqa: D102
        return self._simulation.potential_energy * self._lammps_to_narupa.energy

    @classmethod
    def from_file(
        cls, filename: str, *, units: Optional[UnitSystem] = None, command_line: Optional[List[str]] = None
    ) -> LAMMPSDynamics:
        """
        Load LAMMPS simulation from a file.

        This automatically runs the 'atom_modify map yes' command required to use atom
        IDs and hence allow certain operations such as setting positions.
        :param filename: Filename of LAMMPS input.
        :param units: If the file uses LJ units, a UnitSystem must be provided.
        :return: LAMMPS simulation based on file.
        """
        simulation = LAMMPSSimulation.from_file(filename, units=units, command_line=command_line)
        return cls(simulation)


class LAMMPSIMDFeature(InteractionFeature[LAMMPSDynamics]):
    """Interactive molecular dynamics for use with LAMMPS."""

    def __init__(self, dynamics: LAMMPSDynamics):
        super().__init__(dynamics)

    def calculate_imd(self) -> Tuple[float, Optional[np.ndarray], Optional[np.ndarray]]:
        """Calculate the energy, forces and torques to be applied."""
        if not self.current_interactions:
            return 0, None, None
        energy = 0.0
        forces = np.zeros(shape=(self._system_size, 3))
        torques = np.zeros(shape=(self._system_size, 3))
        has_torques = False
        for interaction in self.current_interactions.values():
            interaction.mark_positions_dirty()
            energy += interaction.potential_energy
            forces[interaction.particle_indices] += interaction.forces
            try:
                torques[interaction.particle_indices] += interaction.torques
                has_torques = True
            except AttributeError:
                pass
        return (
            energy * self.dynamics.simulation._narupa_to_lammps.energy,
            forces * self.dynamics.simulation._narupa_to_lammps.force,
            torques * self.dynamics.simulation._narupa_to_lammps.torque
            if has_torques
            else None,
        )

    @property
    def _system_size(self) -> int:
        return len(self.dynamics.simulation)
