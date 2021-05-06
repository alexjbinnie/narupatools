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

"""Generalization of a server which can broadcast simulations and trajectories."""

from __future__ import annotations

import time
from abc import ABCMeta
from contextlib import contextmanager
from types import TracebackType
from typing import Any, Generator, Generic, Optional, Protocol, Type, TypeVar

from infinite_sets import InfiniteSet
from narupa.app import NarupaImdApplication
from narupa.app.app_server import DEFAULT_NARUPA_PORT
from narupa.core import DEFAULT_SERVE_ADDRESS, NarupaServer
from narupa.essd import DiscoveryServer
from narupa.trajectory import FrameData
from narupa.trajectory.frame_server import (
    PAUSE_COMMAND_KEY,
    PLAY_COMMAND_KEY,
    RESET_COMMAND_KEY,
    STEP_COMMAND_KEY,
)

from narupatools.core.broadcastable import Broadcastable, Broadcaster
from narupatools.core.dynamics import SimulationDynamics
from narupatools.core.event import Event
from narupatools.core.health_check import HealthCheck
from narupatools.core.playable import Playable
from narupatools.frame import KineticEnergy, ParticlePositions, PotentialEnergy
from narupatools.frame.frame_producer import FrameProducer
from narupatools.frame.frame_source import FrameSource
from narupatools.state.view.wrappers import SharedStateServerWrapper
from .shared_state import SharedStateView

TTarget = TypeVar("TTarget")

TTarget_Co = TypeVar("TTarget_Co", contravariant=True)


class OnTargetChanged(Protocol[TTarget_Co]):
    """Callback for when the target of a session is altered."""

    def __call__(
        self, target: Optional[TTarget_Co], previous_target: Optional[TTarget_Co]
    ) -> None:
        """
        Called when the target of a session is altered.

        :param target: The new target of the session, if any.
        :param previous_target: The previous target of the session, if any.
        """
        pass


class Session(Broadcaster, Generic[TTarget], HealthCheck, metaclass=ABCMeta):
    """
    Narupa server bundled with the ability to change what is being displayed.

    A session has one target at a time, which can be a simulation, trajectory or a
    single frame.

    The session has a background loop which periodically gets a FrameData to send to the
    clients. If the target implements FrameSource, this indicates the target can produce
    a frame and this will be sent to the clients. This background loop only sends fields
    which have been marked as dirty. Unlike Narupa, this means the sending of frames is
    completely detached from the MD loop, and hence can be specified in real time.
    """

    def __init__(self, **kwargs: Any):
        """
        Create a session.

        :param kwargs: Parameters to initialise the server with.
        """
        self._target: Optional[TTarget] = None

        self.frame_index = 0
        self._server = self._initialise_server(**kwargs)

        self._frame_producer = FrameProducer(self._produce_frame)
        self._frame_producer.on_frame_produced.add_callback(self._on_frame_produced)
        self._frame_producer.start(block=False)

        self._shared_state = SharedStateView(
            SharedStateServerWrapper(self._server.server)
        )

        narupa_server = self._server.server

        narupa_server._state_service.state_dictionary.content_updated.add_callback(
            self._shared_state.on_dictionary_update
        )

        self._on_target_changed = Event(OnTargetChanged)

        narupa_server.register_command(PLAY_COMMAND_KEY, self.play)  # type: ignore
        narupa_server.register_command(RESET_COMMAND_KEY, self.restart)  # type: ignore
        narupa_server.register_command(STEP_COMMAND_KEY, self.step)  # type: ignore
        narupa_server.register_command(PAUSE_COMMAND_KEY, self.pause)  # type: ignore

    @property
    def app(self) -> NarupaImdApplication:
        """Underlying Narupa application."""
        return self._server

    @property
    def shared_state(self) -> SharedStateView:
        """Shared state of the session."""
        return self._shared_state

    def _on_frame_produced(self, *, frame: FrameData, **kwargs: Any) -> None:
        self._server.frame_publisher.send_frame(self.frame_index, frame)
        self.frame_index += 1

    def run(self, *, block: bool = False) -> None:
        """
        Run the target if it is a playable object.

        :raises RuntimeError: The target is not an instance of Playable.
        :param block: Should this be run in blocking mode.
        """
        if not isinstance(self._target, Playable):
            raise RuntimeError("Runner has invalid target.")
        self._target.run(block)

    def stop(self) -> None:
        """Stop the target if it is a playable object."""
        if isinstance(self._target, Playable):
            self._target.stop(wait=True)

    def restart(self) -> None:
        """Restart the target if it is a playable object."""
        if isinstance(self._target, Playable):
            self._target.restart()

    def step(self) -> None:
        """Step the target if it is a playable object."""
        if isinstance(self._target, Playable):
            self._target.step()

    def play(self) -> None:
        """Play the target if it is a playable object."""
        if isinstance(self._target, Playable):
            self._target.play()

    def pause(self) -> None:
        """Pause the target if it is a playable object."""
        if isinstance(self._target, Playable):
            self._target.pause()

    def show(self, target: TTarget, /) -> None:
        """
        Broadcast an object, replacing the previous target if present.

        :param target: Object to start broadcasting.
        """
        self.target = target

    @property
    def target(self) -> Optional[TTarget]:
        """Target which is being broadcast by this session."""
        return self._target

    @target.setter
    def target(self, value: TTarget) -> None:
        if isinstance(self._target, Playable):
            self._target.stop(wait=True)
        if isinstance(self._target, SimulationDynamics):
            self._target.on_reset.remove_callback(self._on_target_reset)
            self._target.on_post_step.remove_callback(self._on_target_step)
        if isinstance(self._target, Broadcastable):
            self._target.end_broadcast(self)

        previous_target = self._target
        self._target = value

        if isinstance(self._target, SimulationDynamics):
            self._target.on_reset.add_callback(self._on_target_reset)
            self._target.on_post_step.add_callback(self._on_target_step)
        if isinstance(self._target, Broadcastable):
            self._target.start_broadcast(self)

        # Trigger a full frame reset
        self._frame_producer.mark_dirty()

        self._on_target_changed.invoke(
            target=self._target, previous_target=previous_target  # type: ignore
        )

    def _on_target_reset(self, **kwargs: Any) -> None:
        self._frame_producer.mark_dirty()

    def _on_target_step(self, **kwargs: Any) -> None:
        self._frame_producer.mark_dirty(
            {ParticlePositions.key, PotentialEnergy.key, KineticEnergy.key}
        )

    def _produce_frame(self, fields: InfiniteSet[str]) -> FrameData:
        """
        Called by the frame producing loop to request a new frame.

        Also a level 0 cantrip.

        :param fields: The set of fields that the frame is interested in.
        :return: A FrameData to send to the clients. If there isn't a target or the
                 target doesn't produce frames, this returns a blank frame.
        """
        if isinstance(self._target, FrameSource):
            return self._target.get_frame(fields=fields)
        return FrameData()

    def _initialise_server(
        self,
        *,
        name: Optional[str] = None,
        address: Optional[str] = None,
        port: Optional[int] = None,
        run_discovery: bool = True,
        discovery_port: Optional[int] = None,
    ) -> NarupaImdApplication:

        address = address or DEFAULT_SERVE_ADDRESS
        if port is None:
            port = DEFAULT_NARUPA_PORT
        server = NarupaServer(address=address, port=port)
        discovery: Optional[DiscoveryServer] = None
        if run_discovery:
            discovery = DiscoveryServer(broadcast_port=discovery_port)
        return NarupaImdApplication(server, discovery, name)

    def close(self) -> None:
        """Remove the current target and close down the server and background tasks."""
        self.target = None
        self._frame_producer.stop(wait=True)
        self._server.close()

    def health_check(self) -> None:  # noqa: D102
        self._frame_producer.health_check()
        if isinstance(self._target, HealthCheck):
            self._target.health_check()

    @classmethod
    @contextmanager
    def start(cls, **kwargs: Any) -> Generator[Session, None, None]:
        """
        Start a session as a context manager, calling close() after completion.

        :param kwargs: Arguments to pass to the session.
        """
        with cls(**kwargs) as session:
            yield session

    def __enter__(self) -> Session:
        return self

    def __exit__(
        self,
        t: Optional[Type[BaseException]],
        value: Optional[BaseException],
        traceback: Optional[TracebackType],
    ) -> None:
        self.close()

    @property
    def address(self) -> str:
        """URL or IP address the server is broadcasting at."""
        return self._server.address

    @property
    def port(self) -> int:
        """Port the session is broadcasting on."""
        return self._server.port

    @property
    def name(self) -> str:
        """Name of the session."""
        return self._server.name

    def start_loop(self) -> None:
        """
        Run an infinite loop waiting for a keyboard interruption.

        This checks the health of the session every second to catch any exceptions
        occurring on one of the background threads.
        """
        while True:
            try:
                time.sleep(1)
                self.health_check()
            except KeyboardInterrupt:
                return
