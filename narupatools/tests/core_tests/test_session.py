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

import time
from typing import Generator

import numpy as np
import pytest

from narupatools.app import Client, Session
from narupatools.core import Playable
from narupatools.util.timing import wait_for


class TestPlayable(Playable):
    def _restart(self) -> None:
        pass

    def _advance(self) -> bool:
        return True


@pytest.fixture
def session() -> Generator[Session, None, None]:
    with Session(port=0, run_discovery=False) as session:
        yield session


@pytest.fixture
def client(session) -> Generator[Client, None, None]:
    with Client.connect_to_session(session) as client:
        client.subscribe_multiplayer()
        yield client


def test_session_initial_target():
    playable = TestPlayable()
    session = Session(playable, port=0, run_discovery=False)
    assert session.target is playable
    session.close()


def test_session_autoplay():
    playable = TestPlayable()
    assert not playable.is_running
    session = Session(playable, port=0, run_discovery=False)
    wait_for(lambda: playable.is_running)
    session.close()
    playable.stop()


def test_session_no_autoplay():
    playable = TestPlayable()
    assert not playable.is_running
    session = Session(playable, autoplay=False, port=0, run_discovery=False)
    assert not playable.is_running
    assert session.target is playable
    session.close()
    playable.stop()


def test_shared_state_numpy_array(session, client):
    session.shared_state["test"] = np.array([0.0, 1.0])
    time.sleep(0.5)
    assert client.shared_state["test"].snapshot() == pytest.approx([0.0, 1.0])

    client.shared_state["test2"] = np.array([1.0, 0.0])
    time.sleep(0.5)
    assert session.shared_state["test2"].snapshot() == pytest.approx([1.0, 0.0])


def test_shared_state_set(session, client):
    session.shared_state["test"] = {0, 1}
    time.sleep(0.5)
    assert client.shared_state["test"].snapshot() == pytest.approx([0.0, 1.0])

    client.shared_state["test2"] = {0, 1}
    time.sleep(0.5)
    assert session.shared_state["test2"].snapshot() == pytest.approx([0.0, 1.0])
