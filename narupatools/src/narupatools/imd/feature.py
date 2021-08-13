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
from typing import Any, Dict, Generic, List, Mapping, Optional, Set, TypeVar

import numpy as np
from narupa.imd import ParticleInteraction
from typing_extensions import Protocol

from narupatools.core.dynamics import (
    OnPostStepCallback,
    OnPreStepCallback,
    OnResetCallback,
)
from narupatools.core.event import Event, EventListener

from .interaction import Interaction
from .interaction_source import InteractionsSource


class OnStartInteractionCallback(Protocol):
    """Callback for when an interaction is first applied to the system."""

    def __call__(self, *, key: str, interaction: ParticleInteraction) -> None:
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


class DynamicsSupportsInteractions(Protocol):
    """Protocol for objects that support interactions."""

    @property
    def on_pre_step(self) -> EventListener[OnPreStepCallback]:  # noqa: D102
        ...

    @property
    def on_post_step(self) -> EventListener[OnPostStepCallback]:  # noqa: D102
        ...

    @property
    def on_reset(self) -> EventListener[OnResetCallback]:  # noqa: D102
        ...

    @property
    def elapsed_time(self) -> float:  # noqa: D102
        ...


TDynamics = TypeVar("TDynamics", bound=DynamicsSupportsInteractions)
TInteraction = TypeVar("TInteraction", bound=Interaction)


class DynamicInteractions:
    """
    Polls a source of interactions before each simulation step.

    This is a wrapper around an interaction source , which is checked before each step
    to add/remove interactions that have appeared or disappeared since the last check.
    """

    def __init__(self, imd: InteractionFeature, source: InteractionsSource):
        self._imd = imd
        self._source = source
        self._interactions: Set[str] = set()

    def update(self, has_reset: bool) -> None:
        """Update the interactions based on the current state of the source."""
        if has_reset:
            self._interactions.clear()
        interactions = self._source.active_interactions
        for key in list(self._interactions):
            if key not in interactions:
                self._imd.remove_interaction(key)
                self._interactions.remove(key)
        for key in interactions.keys():
            if key not in self._interactions:
                self._imd.add_interaction(interactions[key], key=key)
                self._interactions.add(key)
            else:
                self._imd.update_interaction(key, interactions[key])


class InteractionFeature(Generic[TDynamics, TInteraction]):
    """Interactive Molecular Dynamics manager."""

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

        self._current_interactions: Dict[str, TInteraction] = {}

        self._on_start_interaction: Event[OnStartInteractionCallback] = Event()
        self._on_end_interaction: Event[OnEndInteractionCallback] = Event()

        self._has_reset = False

        self._sources: List[DynamicInteractions] = []

        self._elapsed_work = 0.0

    def add_dynamic_interactions_source(self, source: InteractionsSource) -> None:
        """
        Add a source for interactions which are polled before each dynamics step.

        :param source: Source of IMD interactions.
        """
        self._sources.append(DynamicInteractions(self, source))

    def remove_dynamic_interactions_source(self, source: InteractionsSource) -> None:
        """
        Remove a source for interactions which are polled before each dynamics step.

        :param source: Source of IMD interactions.
        """
        for _source in [s for s in self._sources if s._source is source]:
            self._sources.remove(_source)

    @property
    def dynamics(self) -> TDynamics:
        """Dynamics these IMD interactions are applied to."""
        return self._dynamics

    @property
    def current_interactions(self) -> Mapping[str, TInteraction]:
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

        for source in self._sources:
            source.update(self._has_reset)

        self._has_reset = False

        for interaction in self.current_interactions.values():
            interaction.on_pre_step()

    def _on_post_step(self, **kwargs: Any) -> None:
        for interaction in self.current_interactions.values():
            interaction.on_post_step()

    def add_interaction(
        self, interaction: ParticleInteraction, *, key: Optional[str] = None
    ) -> str:
        """
        Add an interaction to the system.

        :param interaction: Interaction to add.
        :param key: Key to store interaction as. If not provided, automatically
                    generated.
        :return: Key that the interaction was assigned to.
        """
        if key is None:
            key = f"_interaction.{uuid.uuid4()}"
        start_time = self._dynamics.elapsed_time
        constraint = self._make_interaction(
            key=key, interaction=interaction, start_time=start_time
        )
        self._current_interactions[key] = constraint
        self._on_start_interaction.invoke(key=key, interaction=interaction)
        return key

    def update_interaction(self, key: str, interaction: ParticleInteraction) -> None:
        """
        Update an interaction with the given key.

        :param key: Key of the interaction to update.
        :param interaction: ParticleInteraction containing data for the interaction.
        :raises KeyError: Interaction with the given key does not exist.
        """
        self._current_interactions[key].interaction = interaction

    def remove_interaction(self, key: str) -> TInteraction:
        """
        Remove an interaction with the given key.

        :param key: Key of the interaction to remove.
        :raises KeyError: Interaction with the given key does not exist.
        :return: Interaction which was removed.
        """
        constraint = self.current_interactions[key]
        del self._current_interactions[key]

        end_time = self._dynamics.elapsed_time
        duration = end_time - constraint.start_time
        self._elapsed_work += constraint.total_work
        self._on_end_interaction.invoke(
            key=key, work_done=constraint.total_work, duration=duration
        )

        return constraint

    def _on_reset(self, **kwargs: Any) -> None:
        self._has_reset = True

    @abstractmethod
    def _make_interaction(
        self, *, key: str, interaction: ParticleInteraction, start_time: float
    ) -> TInteraction:
        pass
