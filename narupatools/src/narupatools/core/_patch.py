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

"""
Monkey patching to allow protobuf serialization to support more options.

This overrides both the validation of if a dict is serializable and the conversion to a
state update for both the state service and the client. Instead of the protobuf
functions for converting to protobuf Struct's, this uses narupatools's implementation,
which supports a wider range of options such as iterables and NumPy arrays.
"""

from typing import Mapping

from narupa.core import narupa_client
from narupa.protocol.state import StateUpdate
from narupa.state import state_service
from narupa.utilities.change_buffers import DictionaryChange

from narupatools.core.protobuf import dictionary_to_protobuf, update_struct
from narupatools.state.typing import Serializable


def _dictionary_change_to_state_update(change: DictionaryChange) -> StateUpdate:
    changes, removals = change

    update = StateUpdate()
    update_struct(update.changed_keys, changes)
    for key in removals:
        update.changed_keys[key] = None

    return update


def _validate_dict_is_serializable(dictionary: Mapping[str, Serializable]) -> None:
    try:
        dictionary_to_protobuf(dictionary)
    except ValueError as e:
        raise TypeError("Data is not serializable with protobuf.") from e


state_service.dictionary_change_to_state_update = _dictionary_change_to_state_update
state_service.validate_dict_is_serializable = _validate_dict_is_serializable

narupa_client.dictionary_change_to_state_update = _dictionary_change_to_state_update
narupa_client.validate_dict_is_serializable = _validate_dict_is_serializable

__all__ = []  # type: ignore
