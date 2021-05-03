import time
from typing import Generator

import numpy as np
import pytest

from narupatools.app import Client, Session


@pytest.fixture
def session() -> Generator[Session, None, None]:
    with Session(port=0) as session:
        yield session


@pytest.fixture
def client(session) -> Generator[Client, None, None]:
    with Client.connect_to_session(session) as client:
        client.subscribe_multiplayer()
        yield client


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
