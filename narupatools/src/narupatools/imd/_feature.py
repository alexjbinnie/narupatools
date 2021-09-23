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

"""Base manager for handling interactions."""

from __future__ import annotations

import uuid
from abc import abstractmethod
from typing import Any, Callable, Dict, Generic, List, Mapping, Optional, TypeVar, Union

import numpy as np
from typing_extensions import Protocol

from narupatools.core.dynamics import SimulationDynamics
from narupatools.core.event import Event, EventListener
from narupatools.imd.interactions._interaction import Interaction
from narupatools.imd.interactions._parameters import InteractionParameters
from narupatools.physics.thermodynamics import maxwell_boltzmann_velocities


class OnStartInteractionCallback(Protocol):
    """Callback for when an interaction is first applied to the system."""

    def __call__(self, *, key: str, interaction: Interaction) -> None:
        """
        Called when an interaction is first applied to the system.

        :param key: Full shared state key of the interaction.
        :param interaction: Interaction object sent by client.
        """
        ...


class OnEndInteractionCallback(Protocol):
    """Callback for when an interaction is removed from the system."""

    def __call__(self, *, key: str, work_done: float, duration: float) -> None:
        """
        Called when an interaction is removed from the system.

        :param key: Full shared state key of the interaction.
        :param work_done: Total work done by interaction in kilojoules per mole.
        :param duration: Duration of interaction in simulation time in picoseconds.
        """
        ...


TDynamics = TypeVar("TDynamics", bound=SimulationDynamics)


InteractionSource = Union[
    Mapping[str, InteractionParameters],
    Callable[..., Mapping[str, InteractionParameters]],
]


class InteractionFeature(Generic[TDynamics]):
    """
    Interactive Molecular Dynamics manager.

    Depending on the underlying dynamics, the forces applied during IMD need to be applied in different ways.
    This manager takes interactions from one or more sources (as stored as protobuf definitions) and creates
    persistent representations which are subclasses of :class:`Interaction`.
    """

    def __init__(self, dynamics: TDynamics):
        """
        Create a new interaction manager.

        :param dynamics: SimulationDynamics object which notifies the ASEImdProvider
                         when dynamics steps occur.
        """
        self._dynamics = dynamics
        self._dynamics.on_pre_step.add_callback(self._on_pre_step)
        self._dynamics.on_post_step.add_callback(self._on_post_step)
        self._dynamics.on_reset.add_callback(self._on_reset)

        self._current_interactions: Dict[str, Interaction] = {}

        self._on_start_interaction: Event[OnStartInteractionCallback] = Event()
        self._on_end_interaction: Event[OnEndInteractionCallback] = Event()

        self._has_reset = False

        self._sources: List[InteractionSource] = []

        self._user_interaction_keys: List[str] = []

        self._elapsed_work = 0.0

    def add_source(self, source: InteractionSource) -> None:
        """
        Add a source for interactions which are polled before each dynamics step.

        :param source: Source of IMD interactions.
        """
        self._sources.append(source)

    def remove_source(self, source: InteractionSource) -> None:
        """
        Remove a source for interactions which are polled before each dynamics step.

        :param source: Source of IMD interactions.
        """
        for _source in [s for s in self._sources if s is source]:
            self._sources.remove(_source)

    def _source_interactions(self) -> Mapping[str, InteractionParameters]:
        """Collate all interaction sources into one unified dictionary."""
        dict_: Dict[str, InteractionParameters] = {}
        for source in self._sources:
            if callable(source):
                dict_.update(source())
            else:
                dict_.update(source)
        return dict_

    @property
    def dynamics(self) -> TDynamics:
        """Dynamics these IMD interactions are applied to."""
        return self._dynamics

    @property
    def current_interactions(self) -> Mapping[str, Interaction]:
        """Key-index set of interaction constraints currently applied to the system."""
        return self._current_interactions

    @property
    def on_start_interaction(self) -> EventListener[OnStartInteractionCallback]:
        """Event triggered when an interaction is first applied to the system."""
        return self._on_start_interaction

    @property
    def on_end_interaction(self) -> EventListener[OnEndInteractionCallback]:
        """Event triggered when an interaction is removed from the system."""
        return self._on_end_interaction

    @property
    def total_work(self) -> float:
        """
        Total work performed by all active interactions, in kilojoules per mole.

        This does not include interactions which have now finished.
        """
        total_work = [
            interaction.total_work for interaction in self.current_interactions.values()
        ]
        return self._elapsed_work + np.array(total_work).sum()  # type:ignore

    @property
    def work_last_step(self) -> float:
        """
        Total work performed last step by all active interactions.

        The returned work is in kilojoules per mole.
        """
        work_last_step = [
            interaction.work_last_step
            for interaction in self.current_interactions.values()
        ]
        return np.array(work_last_step).sum()  # type:ignore

    @property
    def potential_energy(self) -> float:
        """
        Total potential energy from all active interactions.

        The returned work is in kilojoules per mole.
        """
        potential_energy = [
            interaction.potential_energy
            for interaction in self.current_interactions.values()
        ]
        return np.array(potential_energy).sum()  # type:ignore

    @property
    def imd_forces(self) -> np.ndarray:
        """
        Total interactive forces from all active interactions.

        This returns a (N, 3) NumPy array, where N is the number of atoms in the system.
        The forces are in kilojoules per mole per nanometer.
        """
        forces = np.zeros((self._system_size, 3))
        for interaction in self._current_interactions.values():
            forces[interaction.particle_indices, :] += interaction.forces
        return forces

    @property
    @abstractmethod
    def _system_size(self) -> int:
        pass

    def close(self) -> None:
        """Remove this ASEImdProvider from the dynamics it was applied to."""
        self._dynamics.on_pre_step.remove_callback(self._on_pre_step)
        self._dynamics.on_post_step.remove_callback(self._on_post_step)
        self._dynamics.on_reset.remove_callback(self._on_reset)

    def _on_pre_step(self, **kwargs: Any) -> None:
        if self._has_reset:
            for key in list(self.current_interactions.keys()):
                self.remove_interaction(key)
            self._user_interaction_keys.clear()

        _source_interactions = self._source_interactions()
        for key in list(self.current_interactions.keys()):
            if (
                key not in _source_interactions
                and key not in self._user_interaction_keys
            ):
                self.remove_interaction(key)
        for key, interaction in _source_interactions.items():
            if key in self._user_interaction_keys:
                continue
            if key not in self.current_interactions:
                self._add_interaction(key, interaction)
            else:
                self.update_interaction(key=key, interaction=interaction)

        self._has_reset = False

        for interaction in self.current_interactions.values():
            interaction.on_pre_step()

    def _on_post_step(self, **kwargs: Any) -> None:
        for interaction in self.current_interactions.values():
            interaction.on_post_step(**kwargs)

    def add_interaction(
        self,
        /,
        interaction_data: InteractionParameters,
        *,
        key: Optional[str] = None,
    ) -> str:
        """
        Add an interaction to the system.

        :param interaction_data: Interaction to add.
        :param key: Key to store interaction as. If not provided, automatically
                    generated.
        :return: Key that the interaction was assigned to.
        """
        if key is None:
            key = f"_interaction.{uuid.uuid4()}"
        self._user_interaction_keys.append(key)
        self._add_interaction(key, interaction_data)
        return key

    def _add_interaction(
        self,
        key: str,
        interaction_data: InteractionParameters,
    ) -> None:
        """
        Add an interaction to the system.

        :param interaction_data: Interaction to add.
        :param key: Key to store interaction as. If not provided, automatically
                    generated.
        :return: Key that the interaction was assigned to.
        """
        start_time = self._dynamics.elapsed_time
        interaction = self.create_interaction(
            key=key, start_time=start_time, interaction=interaction_data
        )
        self._current_interactions[key] = interaction
        self._on_start_interaction.invoke(key=key, interaction=interaction)

    def update_interaction(
        self, *, key: str, interaction: InteractionParameters
    ) -> None:
        """
        Update an interaction with the given key.

        :param key: Key of the interaction to update.
        :param interaction: ParticleInteraction containing data for the interaction.
        :raises KeyError: Interaction with the given key does not exist.
        """
        self._current_interactions[key].update(interaction)

    def create_interaction(
        self, *, key: str, interaction: InteractionParameters, start_time: float
    ) -> Interaction:
        """
        Create a new interaction from interaction data.

        This allows subclasses to override how interactions are created. At the very
        least, it is strongly encouraged to still use :meth:`Interaction.create`, as that
        creates the correct subclass depending on the interaction type described in the
        provided :class:`InteractionData`.

        :param key: Key of the interaction.
        :param interaction: Interaction data to base interaction on.
        :param start_time: Start time of the interaction in picoseconds.
        :return:
        """
        return Interaction.create(
            key=key,
            interaction=interaction,
            start_time=start_time,
            dynamics=self.dynamics,
        )

    def clear_interactions(self) -> None:
        """Remove all current interactions."""
        keys = self.current_interactions.keys()
        for key in list(keys):
            self.remove_interaction(key)

    def remove_interaction(self, key: str) -> Interaction:
        """
        Remove an interaction with the given key.

        :param key: Key of the interaction to remove.
        :raises KeyError: Interaction with the given key does not exist.
        :return: Interaction which was removed.
        """
        interaction = self.current_interactions[key]
        del self._current_interactions[key]

        end_time = self._dynamics.elapsed_time
        duration = end_time - interaction.start_time
        self._elapsed_work += interaction.total_work
        self._on_end_interaction.invoke(
            key=key, work_done=interaction.total_work, duration=duration
        )

        if key in self._user_interaction_keys:
            self._user_interaction_keys.remove(key)

        if interaction.velocity_reset:
            velocities = self.dynamics.velocities
            velocities[interaction.particle_indices] = maxwell_boltzmann_velocities(
                masses=self.dynamics.masses[interaction.particle_indices],
                temperature=300,
            )
            self.dynamics.velocities = velocities  # type: ignore

        return interaction

    def _on_reset(self, **kwargs: Any) -> None:
        self._has_reset = True
