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
from contextlib import contextmanager

from narupa.app import NarupaImdClient
from narupa.core import NarupaServer
from narupa.utilities.change_buffers import DictionaryChange

from narupatools.state.view._wrappers import (
    SharedStateClientWrapper,
    SharedStateServerWrapper,
)


class SharedStateClientWrapperWithDelay(SharedStateClientWrapper):
    """
    Shared state view with delay, to ensure changes have been applied before testing
    """

    def __setitem__(self, key, value):
        super().__setitem__(key, value)
        time.sleep(0.1)

    def __delitem__(self, key):
        super().__delitem__(key)
        time.sleep(0.1)


@contextmanager
def dictionary_shared_state_wrapper(dictionary):
    yield dictionary


@contextmanager
def narupa_server_with_initial_state(dictionary):
    server = NarupaServer(address="localhost", port=0)
    server.update_state(None, DictionaryChange(updates=dictionary))
    yield server
    server.close()


@contextmanager
def narupa_client_with_initial_state(dictionary):
    server = NarupaServer(address="localhost", port=0)
    client = NarupaImdClient(multiplayer_address=(server.address, server.port))
    client.subscribe_multiplayer()
    server.update_state(None, DictionaryChange(updates=dictionary))
    time.sleep(0.1)
    yield client
    client.close()
    server.close()


@contextmanager
def server_shared_state_wrapper(dictionary):
    with narupa_server_with_initial_state(dictionary) as server:
        yield SharedStateServerWrapper(server)


@contextmanager
def client_shared_state_wrapper(dictionary):
    with narupa_client_with_initial_state(dictionary) as client:
        yield SharedStateClientWrapperWithDelay(client)


# Possible views of a shared state dictionary, namely a raw dictionary, a server and a
# client
SHARED_STATE_WRAPPERS = [
    dictionary_shared_state_wrapper,
    server_shared_state_wrapper,
    client_shared_state_wrapper,
]
