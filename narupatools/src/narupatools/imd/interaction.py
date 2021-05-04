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

"""Base implementation of an interaction."""

from abc import ABCMeta, abstractmethod
from typing import Generic, Optional, Tuple, TypeVar

import numpy as np
from narupa.imd import ParticleInteraction

from narupatools.physics.typing import Vector3Array

_TDynamics = TypeVar("_TDynamics")


class Interaction(Generic[_TDynamics], metaclass=ABCMeta):
    """Base class for a stateful interaction that tracks work done."""

    def __init__(
        self,
        *,
        dynamics: _TDynamics,
        key: str,
        interaction: ParticleInteraction,
        start_time: float,
    ):
        """
        Create a new interaction.

        :param dynamics: Underlying dynamics this interaction is affecting.
        :param key: Unique key to identify this interaction.
        :param interaction: Initial parameters of the interaction.
        :param start_time: Start time of interaction in simulation time in picoseconds
        """
        self.key: str = key
        self._dynamics: _TDynamics = dynamics
        self._total_work: float = 0.0
        self._last_step_work: Optional[float] = None
        self._interaction: ParticleInteraction = interaction
        self._previous_positions: Vector3Array = np.zeros(0)
        self._previous_forces: Vector3Array = np.zeros(0)
        self._start_time: float = start_time

        self._energy: Optional[float] = None
        self._forces: Optional[Vector3Array] = None

        self._c = 0.0

    @property
    def particle_indices(self) -> np.ndarray:
        """List of indices affected by this interaction."""
        return self._interaction.particles

    @property
    def start_time(self) -> float:
        """Start time of the interaction in picoseconds."""
        return self._start_time

    @property
    def total_work(self) -> float:
        """Total work performed by interaction in kilojoules per mole."""
        return self._total_work

    @property
    def work_last_step(self) -> float:
        """Work performed last step in kilojoules per mole."""
        if self._last_step_work is None:
            raise AttributeError("No dynamics steps have occurred yet.")
        return self._last_step_work

    @property
    def interaction(self) -> ParticleInteraction:
        """Current parameters of the interaction."""
        return self._interaction

    @interaction.setter
    def interaction(self, value: ParticleInteraction) -> None:
        self._interaction = value

    @property
    def potential_energy(self) -> float:
        """Potential energy of the interaction, in kilojoules per mole."""
        if self._energy is None:
            self._forces, self._energy = self.update_energy_and_forces()
        return self._energy

    @property
    def forces(self) -> np.ndarray:
        """
        Forces that will be applied by the interaction.

        The forces are in kilojoules per mole per nanometer.

        This is a (N, 3) NumPy array, where N is the number of particles affected by
        this interaction.
        """
        if self._forces is None:
            self._forces, self._energy = self.update_energy_and_forces()
        return self._forces

    def on_pre_step(self) -> None:
        """Perform any tasks necessary before a dynamics step."""
        self._previous_positions = self.get_positions()
        self._forces, self._energy = self.update_energy_and_forces()
        self._previous_forces = self.forces

    def on_post_step(self) -> None:
        """Perform any tasks necessary after a dynamics step."""
        _current_positions = self.get_positions()

        self._forces, self._energy = self.update_energy_and_forces()

        work_this_step = 0.0

        for i in range(len(self._interaction.particles)):
            # Use trapezoidal rule to calculate single step of integral F.dS
            F = 0.5 * (self._previous_forces[i] + self.forces[i])
            dS = _current_positions[i] - self._previous_positions[i]
            work_this_step += np.dot(F, dS)

        self._last_step_work = work_this_step
        self._total_work += work_this_step

        self._previous_positions = _current_positions

    @abstractmethod
    def update_energy_and_forces(self) -> Tuple[Vector3Array, float]:
        """
        Calculate the forces and energy of this interaction.

        Overriding this allows a subclass to implement features such as caching.

        :return: Tuple of forces (in kilojoules per mole per nanometer) and energy (in
                 kilojoules per mole).
        """
        pass

    @abstractmethod
    def get_positions(self) -> Vector3Array:
        """
        Get the positions of the affected particles for this interaction, in nanometers.

        These are needed to calculate work.

        :return: Array of positions for the particles affected by this interaction.
        """
        pass
