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

"""Base class for simulations that support IMD."""

from abc import ABCMeta, abstractmethod
from typing import Any

from narupatools.app import Broadcastable, Session
from narupatools.core.dynamics import SimulationDynamics
from narupatools.imd._feature import InteractionFeature


class InteractiveSimulationDynamics(
    SimulationDynamics, Broadcastable, metaclass=ABCMeta
):
    """
    Base class for any simulation dynamics which supports IMD.

    Dynamics that support IMD perform the following additional features:

    * Allow interactions to be applied via an :obj:`InteractionFeature` interface.
    * If added to a :obj:`Session`, reads interactions from the shared state. Interactions
      consist of any keys starting with 'interaction.'
    * Add interaction feedback items to the shared state, based on the current interaction.
      This feedback contains realtime information on the interaction such as total work.
    """

    @property
    @abstractmethod
    def imd(self) -> InteractionFeature:
        """Access to the interactions currently being applied to the simulation."""
        raise NotImplementedError

    def _interaction_ended(self, *, key: str, **kwargs: Any) -> None:
        """Clear up interaction feedback based on that interaction."""
        if key.startswith("interaction."):
            feedback_key = "interaction_feedback." + key[12:]
            del self._shared_state[feedback_key]

    def _send_interaction_feedback(self, **kwargs: Any) -> None:
        for key, interaction in self.imd.current_interactions.items():
            feedback_key = "interaction_feedback." + key[12:]
            feedback = interaction.create_feeback()
            self._shared_state[feedback_key] = feedback

    def start_broadcast(self, session: Session) -> None:  # noqa: D102
        self._shared_state = session.shared_state
        self.imd.add_source(session.shared_state.interactions.snapshot)
        self.imd.on_end_interaction.add_callback(self._interaction_ended)
        # Low priority to ensure interaction has updated first
        self.on_post_step.add_callback(self._send_interaction_feedback, priority=-10)

    def end_broadcast(self, session: Session) -> None:  # noqa: D102
        self.imd.remove_source(session.shared_state.interactions.snapshot)
        self.imd.on_end_interaction.remove_callback(self._interaction_ended)
