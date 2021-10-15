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

from typing import Any, Dict, Optional

import numpy as np
from infinite_sets import InfiniteSet, everything
from narupa.trajectory import FrameData

from narupatools.imd import InteractiveSimulationDynamics, SetAndClearInteractionFeature
from narupatools.physics.typing import (
    ScalarArray,
    Vector3,
    Vector3Array,
    Vector3ArrayLike,
)
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

        simulation.add_imd_force()

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

    def _get_frame(self, fields: InfiniteSet[str]) -> FrameData:
        return self._simulation.get_frame(fields=fields)

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
        cls, filename: str, units: Optional[UnitSystem] = None
    ) -> LAMMPSDynamics:
        """
        Load LAMMPS simulation from a file.

        This automatically runs the 'atom_modify map yes' command required to use atom
        IDs and hence allow certain operations such as setting positions.
        :param filename: Filename of LAMMPS input.
        :param units: If the file uses LJ units, a UnitSystem must be provided.
        :return: LAMMPS simulation based on file.
        """
        simulation = LAMMPSSimulation.from_file(filename, units)
        return cls(simulation)


class LAMMPSIMDFeature(SetAndClearInteractionFeature[LAMMPSDynamics]):
    """Interactive molecular dynamics for use with LAMMPS."""

    def __init__(self, dynamics: LAMMPSDynamics):
        super().__init__(dynamics)

    def _on_pre_step(self, **kwargs: Any) -> None:
        super()._on_pre_step(**kwargs)
        self._calculate_and_apply_interactions()

    def _set_forces(self, forces: Dict[int, Vector3], /) -> None:
        for index, force in forces.items():
            self._dynamics.simulation.set_imd_force(
                index, force * self.dynamics._lammps_to_narupa.force
            )

    def _clear_forces(self, indices: InfiniteSet[int] = everything(), /) -> None:
        for index in set(range(self._system_size)) & indices:  # type: ignore[operator]
            self._dynamics.simulation.clear_imd_force(index)

    @property
    def _system_size(self) -> int:
        return len(self.dynamics.simulation)
