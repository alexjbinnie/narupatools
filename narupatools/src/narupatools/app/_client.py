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

from infinite_sets import InfiniteSet
from narupa.app import NarupaImdClient
from narupa.app.app_server import DEFAULT_NARUPA_PORT
from narupa.core import DEFAULT_CONNECT_ADDRESS
from narupa.imd import ParticleInteraction
from narupa.imd.imd_state import dict_to_interaction
from narupa.trajectory import FrameData

from narupatools.core.event import Event, EventListener
from narupatools.imd.interactions import InteractionParameters
from narupatools.override import override
from narupatools.state.view import SharedStateClientWrapper

from ._session import Session
from ._shared_state import SessionSharedState, SharedStateMixin
from ..frame import FrameSource


class OnFrameReceivedCallback(Protocol):
    """Callback for when a frame is received by a client."""

    def __call__(self, *, changes: AbstractSet[str]) -> None:  # noqa: D102
        ...


class Client(NarupaImdClient, SharedStateMixin):
    """
    An extension of the standard Narupa python client with more features.

    These additional features include:

    * A callback :obj:`on_frame_received` which is invoked each time a new frame is received.
    * A :obj:`shared_state` field which exposes the shared state through the narupatools shared state API.
    * Allows interactions to be started using a narupatools :obj:`InteractionData` object in addition to
      the :obj:`ParticleInteraction` object.
    * Fixes a bug where the client starts receiving frames before it has finished setting up the necessary
      variables to store them.
    """

    @override(FrameSource.get_frame)
    def get_frame(self, *, fields: InfiniteSet[str]) -> FrameData:  # noqa: D102
        return self.current_frame  # type: ignore

    _on_frame_received_event: Event[OnFrameReceivedCallback]

    def __init__(self, **kwargs: Any) -> None:
        """
        Create a new :obj:`Client`.

        It is recommended to use the class methods :obj:`connect_to_session` and
        :obj:`connect_to_server` to create a new client, rather than using this constructor directly::

            with Client.connect_to_session(session) as client:
                ...

            with Client.connect_to_server(address="localhost", port=38801) as client:
                ...

        :param kwargs: Keyword arguments to pass to :obj:`NarupaImdClient`.
        """
        self._on_frame_received_event = Event()

        self._shared_state = SessionSharedState(SharedStateClientWrapper(self))

        super().__init__(**kwargs)

        self._multiplayer_client._state.content_updated.add_callback(  # type: ignore
            self._shared_state._on_dictionary_update
        )

    @property
    def shared_state(self) -> SessionSharedState:
        """Shared state as currently seen by the client."""
        return self._shared_state

    @property
    def on_frame_received(self) -> EventListener[OnFrameReceivedCallback]:
        """Event triggered when a frame is received by the client."""
        return self._on_frame_received_event

    @override(NarupaImdClient._on_frame_received)
    def _on_frame_received(self, frame_index: int, frame: FrameData) -> None:
        changes = frame.array_keys | frame.value_keys
        super()._on_frame_received(frame_index, frame)
        self._on_frame_received_event.invoke(changes=changes)

    @classmethod
    @contextmanager
    def connect_to_session(cls, session: Session, /) -> Generator[Client, None, None]:
        """
        Create a new client and connect to the provided session.

        This is for when the session and the client are both running on the same computer. The
        connection will still go through the gRPC channels, as this method purely uses the
        provided session to get its address and port.

        :param session: Session that is running.
        :yields: New :obj:`Client` which is connected to the given session.
        """
        with cls.connect_to_server(address="localhost", port=session.port) as client:
            yield client

    @classmethod
    @contextmanager
    def connect_to_server(
        cls, *, address: Optional[str] = None, port: Optional[int] = None
    ) -> Generator[Client, None, None]:
        """
        Connect to a server with the given address and port.

        :param address: Address to connect to, or None to use the default.
        :param port: Port to connect to, or None to use the default.
        :yield: New :obj:`Client` connected to all available services on the server at the given
                 destination.
        """
        address = address or DEFAULT_CONNECT_ADDRESS
        port = port or DEFAULT_NARUPA_PORT
        url = (address, port)
        with cls(
            trajectory_address=url, imd_address=url, multiplayer_address=url
        ) as client:
            yield client

    @override(NarupaImdClient.start_interaction)
    def start_interaction(
        self,
        interaction: Optional[Union[InteractionParameters, ParticleInteraction]] = None,
    ) -> str:
        """
        Start an interaction with the IMD server.

        This method can take either a :obj:`ParticleInteraction` (narupa representation of an interaction)
        or :obj:`InteractionParameters` (narupatools representation), with :obj:`InteractionData` being designed
        to support a wider range of interaction types.

        :param interaction: Initial interaction data.
        :return: Unique interaction ID generated for this interaction.
        """
        if isinstance(interaction, InteractionParameters):
            interaction = dict_to_interaction(interaction.serialize())
        return super().start_interaction(interaction)  # type: ignore[no-any-return]
