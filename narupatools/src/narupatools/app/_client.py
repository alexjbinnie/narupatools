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

"""Python client for connecting to Narupa servers."""

from __future__ import annotations

from contextlib import contextmanager
from typing import AbstractSet, Any, Generator, Optional, Protocol, Union

from narupa.app import NarupaImdClient
from narupa.app.app_server import DEFAULT_NARUPA_PORT
from narupa.core import DEFAULT_CONNECT_ADDRESS
from narupa.imd import ParticleInteraction
from narupa.imd.imd_state import dict_to_interaction
from narupa.trajectory import FrameData

from narupatools.core.event import Event, EventListener
from narupatools.imd.interactions._interactiondata import InteractionData
from narupatools.state.view._wrappers import SharedStateClientWrapper

from ._session import Session
from ._shared_state import SessionSharedState


class OnFrameReceivedCallback(Protocol):
    """Callback for when a frame is received by a client."""

    def __call__(self, *, changes: AbstractSet[str]) -> None:  # noqa: D102
        ...


class Client(NarupaImdClient):
    """
    An extension of NarupaImdClient with more features.

    These features include:

    * Callback for when a frame is received.
    * Exposing the shared state through the narupatools shared state API.
    * Fixing a bug where the client connects to frames before it has defined the
      needed variables to store them.

    """

    _on_frame_received_event: Event[OnFrameReceivedCallback]

    def __init__(self, **kwargs: Any) -> None:
        """
        Create a client.

        :param kwargs: Keyword arguments to pass to NarupaImdClient.
        """
        self._on_frame_received_event = Event()

        self._shared_state = SessionSharedState(SharedStateClientWrapper(self))

        super().__init__(**kwargs)

        self._multiplayer_client._state.content_updated.add_callback(  # type: ignore
            self._shared_state.on_dictionary_update
        )

    @property
    def shared_state(self) -> SessionSharedState:
        """Shared state of the client."""
        return self._shared_state

    @property
    def on_frame_received(self) -> EventListener[OnFrameReceivedCallback]:
        """Event triggered when a frame is received by the client."""
        return self._on_frame_received_event

    def _on_frame_received(self, frame_index: int, frame: FrameData) -> None:
        changes = frame.array_keys | frame.value_keys
        super()._on_frame_received(frame_index, frame)
        self._on_frame_received_event.invoke(changes=changes)

    @classmethod
    @contextmanager
    def connect_to_session(cls, session: Session, /) -> Generator[Client, None, None]:
        """
        Connect to a session.

        :param session: Session that is running.
        :yields: Client connected to the given session.
        """
        with cls.connect_to_server(address="localhost", port=session.port) as client:
            yield client

    @classmethod
    @contextmanager
    def connect_to_server(
        cls, *, address: Optional[str] = None, port: Optional[int] = None
    ) -> Generator[Client, None, None]:
        """
        Connect to a server.

        :param address: Address to connect to, or None to use the default.
        :param port: Port to connect to, or None to use the default.
        :yield: Client connected to all available services on the server at the given
                 destination.
        """
        address = address or DEFAULT_CONNECT_ADDRESS
        port = port or DEFAULT_NARUPA_PORT
        url = (address, port)
        with cls(
            trajectory_address=url, imd_address=url, multiplayer_address=url
        ) as client:
            yield client

    def start_interaction(
        self, interaction: Optional[Union[InteractionData, ParticleInteraction]] = None
    ) -> str:
        """
        Start an interaction with the IMD server.

        This method can take either ParticleInteraction (narupa representation of an interaction)
        or InteractionData (narupatools representation).

        :param interaction: An optional :class: ParticleInteraction with which
            to begin.
        :return: The unique interaction ID of this interaction, which can be
            used to update the interaction with
            :func:`~NarupaClient.update_interaction`.
        """
        if isinstance(interaction, InteractionData):
            interaction = dict_to_interaction(interaction.serialize())
        return super().start_interaction(interaction)  # type: ignore[no-any-return]
