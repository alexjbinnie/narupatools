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

from typing import Generator

import pytest
from narupa.app import NarupaImdClient
from narupa.core import NarupaServer

from fixtures import SharedStateClientWrapperWithDelay


class Unserializable:
    pass


KEY = "abc"
VALUE = "value"


@pytest.fixture
def server() -> Generator[NarupaServer, None, None]:
    server = NarupaServer(address="localhost", port=0)
    yield server
    server.close()


@pytest.fixture
def client(server) -> Generator[NarupaImdClient, None, None]:
    client = NarupaImdClient(multiplayer_address=(server.address, server.port))
    client.subscribe_multiplayer()
    yield client
    client.close()


@pytest.fixture
def wrapper(client):
    return SharedStateClientWrapperWithDelay(client)


def test_setitem(server, client, wrapper):
    wrapper[KEY] = VALUE
    assert client.get_shared_value(KEY) == VALUE
    assert server.copy_state()[KEY] == VALUE


def test_getitem(server, client, wrapper):
    wrapper[KEY] = VALUE
    assert wrapper[KEY] == VALUE


def test_delitem(server, client, wrapper):
    wrapper[KEY] = VALUE
    del wrapper[KEY]
    assert KEY not in wrapper


def test_setitem_unserializable(server, client, wrapper):
    with pytest.raises(TypeError):
        wrapper[KEY] = Unserializable()
