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

import pytest

from narupatools.app import Client, Session
from narupatools.core.timing import wait_for
from narupatools.frame import ParticlePositions


@pytest.fixture
def session():
    with Session.start(port=0) as session:
        yield session


@pytest.fixture
def session_villin_openmm(session, villin_ase_simulation_dynamics):
    session.show(villin_ase_simulation_dynamics)
    villin_ase_simulation_dynamics.run(block=False)
    time.sleep(0.5)
    yield session
    villin_ase_simulation_dynamics.stop(wait=True)


@pytest.fixture
def client(session):
    with Client.connect_to_session(session) as client:
        yield client


def test_client_gets_frame(session_villin_openmm, client):
    client.subscribe_to_frames()
    wait_for(lambda: ParticlePositions.key in client.current_frame)
    assert ParticlePositions.key in client.current_frame
    assert ParticlePositions.key in client.current_frame.copy()
