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

from __future__ import annotations

from abc import ABCMeta, abstractmethod
from typing import Any, Dict, Generic, Optional, Type, TypeVar

import numpy as np

from narupatools.physics.typing import Vector3Array

from ._feedback import InteractionFeedback
from ._parameters import InteractionParameters
from ...frame import DynamicStructureProperties

_TInteractionData = TypeVar("_TInteractionData", bound=InteractionParameters)


class Interaction(Generic[_TInteractionData], metaclass=ABCMeta):
    """Base class for a stateful interaction that tracks work done."""

    _types: Dict[str, Type[Interaction]] = {}

    _feedback_type: Type[InteractionFeedback] = InteractionFeedback

    def __init__(
        self,
        *,
        dynamics: DynamicStructureProperties,
        key: str,
        start_time: float,
        interaction: _TInteractionData,
    ):
        """
        Create a new interaction.

        :param dynamics: Underlying dynamics this interaction is affecting.
        :param key: Unique key to identify this interaction.
        :param start_time: Start time of interaction in simulation time in picoseconds
        """
        self.interaction_type: str = interaction.interaction_type
        self.key: str = key
        self._dynamics = dynamics
        self._total_work: float = 0.0
        self._last_step_work: Optional[float] = None
        self._previous_positions: Vector3Array = np.zeros(0)
        self._previous_forces: Vector3Array = np.zeros(0)
        self._start_time: float = start_time

        self.velocity_reset = True
        self.mass_weighted = True

        self._energy: Optional[float] = None
        self._forces: Optional[Vector3Array] = None
        self._torques: Optional[Vector3Array] = None

        self._particles = interaction.particles
        self.update(interaction)

        self._forces_energy_dirty = True

    @classmethod
    def register_interaction_type(
        cls, interaction_type: str, python_type: Type[Interaction], /
    ) -> None:
        """Register a new type of interaction."""
        cls._types[interaction_type] = python_type

    @classmethod
    def _get_class(cls, interaction_type: str) -> Type[Interaction]:
        return cls._types[interaction_type]

    @classmethod
    def create(
        cls,
        *,
        key: str,
        dynamics: DynamicStructureProperties,
        start_time: float,
        interaction: InteractionParameters,
    ) -> Interaction:
        """
        Create a new interaction instance.

        :param key: Key of the interaction.
        :param dynamics: Dynamics used by the interaction.
        :param start_time: Start time of the interaction in picoseconds.
        :param interaction: Initial interaction parameters.
        :return: Interaction instance.
        """
        return cls._get_class(interaction.interaction_type)(
            key=key, dynamics=dynamics, start_time=start_time, interaction=interaction
        )

    def update(self, interaction: _TInteractionData) -> None:
        """Update the interaction based on new data."""
        try:
            self.velocity_reset = interaction.reset_velocities
        except AttributeError:
            self.velocity_reset = False

        try:
            self.mass_weighted = interaction.mass_weighted
        except AttributeError:
            self.mass_weighted = True

    @property
    def particle_indices(self) -> np.ndarray:
        """List of indices affected by this interaction."""
        return self._particles

    @particle_indices.setter
    def particle_indices(self, value: np.ndarray) -> None:
        self._particles = value

    def __len__(self) -> int:
        return len(self._particles)

    @property
    def dynamics(self) -> DynamicStructureProperties:
        """Underlying dynamics this interaction is affecting."""
        return self._dynamics

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
    def potential_energy(self) -> float:
        """Potential energy of the interaction, in kilojoules per mole."""
        if self._forces_energy_dirty:
            self.calculate_forces_and_energy()
            self._forces_energy_dirty = False
        if self._energy is None:
            raise AttributeError
        return self._energy

    @property
    def forces(self) -> np.ndarray:
        """
        Forces that will be applied by the interaction.

        The forces are in kilojoules per mole per nanometer.

        This is a (N, 3) NumPy array, where N is the number of particles affected by
        this interaction.
        """
        if self._forces_energy_dirty:
            self.calculate_forces_and_energy()
            self._forces_energy_dirty = False
        if self._forces is None:
            raise AttributeError
        return self._forces

    @property
    def torques(self) -> np.ndarray:
        """
        Torques that will be applied by the interaction.

        The torques are in kilojoules per mole.

        This is a (N, 3) NumPy array, where N is the number of particles affected by
        this interaction.
        """
        if self._forces_energy_dirty:
            self.calculate_forces_and_energy()
            self._forces_energy_dirty = False
        if self._torques is None:
            raise AttributeError
        return self._torques

    def on_pre_step(self) -> None:
        """Perform any tasks necessary before a dynamics step."""
        self._previous_positions = self._dynamics.positions[self._particles]
        self._previous_forces = self.forces

    def on_post_step(self, timestep: float, **kwargs: Any) -> None:
        """Perform any tasks necessary after a dynamics step."""
        _current_positions = self._dynamics.positions[self._particles]

        new_forces = self.forces

        work_this_step = (
            0.5
            * (
                (self._previous_forces + new_forces)
                * (_current_positions - self._previous_positions)
            ).sum()
        )

        self._last_step_work = work_this_step
        self._total_work += work_this_step

        self._previous_positions = _current_positions

    def create_feeback(self) -> InteractionFeedback:
        """Create a feedback object describing the state of the interaction."""
        feedback = self._feedback_type()
        feedback.interaction_type = self.interaction_type
        feedback.total_work = self._total_work
        feedback.potential_energy = self.potential_energy
        return feedback

    @abstractmethod
    def calculate_forces_and_energy(self) -> None:
        """
        Calculate the forces, torques and energy of this interaction.

        Overriding this allows a subclass to implement features such as caching.
        """

    def mark_positions_dirty(self) -> None:
        """Mark positions as having changed."""
        self._forces_energy_dirty = True

    def mark_velocities_dirty(self) -> None:
        """Mark velocities as having changed."""
        self._forces_energy_dirty = True
