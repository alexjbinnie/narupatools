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

"""Protobuf utilities which support more options than the base Narupa converter."""

from typing import Any, Iterable, Mapping

import six
from google.protobuf.struct_pb2 import ListValue, Struct

from narupatools.physics import quaternion
from narupatools.state import SerializableObject
from narupatools.state.typing import Serializable


class ProtobufSerialisationError(ValueError):
    """Error raised when something cannot be serialized to protobuf."""

    def __init__(self, obj: Any):
        super().__init__(f"Can't convert {obj} to a protobuf representation.")


def dictionary_to_protobuf(mapping: Mapping[str, Serializable]) -> Struct:
    """
    Convert an arbitrary mapping to a protobuf Struct.

    Improved version of :func:`~narupa.utilities.protobuf_utilities.dict_to_struct`,
    that can serialize more complex types such as iterables, mappings and NumPy arrays.

    :param mapping: Dictionary to convert.
    :return: :class:`~google.protobuf.struct_pb2.Struct` containing copies of all the
             items of the dictionary.
    """
    struct = Struct()
    update_struct(struct, mapping)
    return struct


def update_struct(struct: Struct, mapping: Mapping[str, Serializable]) -> None:
    """
    Update a protobuf struct in place from a dictionary-like object.

    :param struct: Struct to update.
    :param mapping: Mapping to provide values to insert into struct.
    """
    for key, value in mapping.items():
        _set_protobuf_value(struct.fields[key], value)


def _extend_listvalue(list_value: ListValue, items: Iterable[Serializable]) -> None:
    for item in items:
        _set_protobuf_value(list_value.values.add(), item)


def _set_protobuf_value(struct_value: Any, value: Serializable) -> None:
    if value is None:
        struct_value.null_value = 0
    elif isinstance(value, bool):
        struct_value.bool_value = value
    elif isinstance(value, quaternion):  # type: ignore
        struct_value.list_value.Clear()
        _extend_listvalue(struct_value.list_value, value.components)
    elif isinstance(value, six.string_types):
        struct_value.string_value = value
    elif isinstance(value, (int, float)):
        struct_value.number_value = value
    elif isinstance(value, Mapping):
        struct_value.struct_value.Clear()
        update_struct(struct_value.struct_value, value)
    elif isinstance(value, (Iterable, ListValue)):
        struct_value.list_value.Clear()
        _extend_listvalue(struct_value.list_value, value)
    elif isinstance(value, SerializableObject):
        _set_protobuf_value(struct_value, value.serialize())
    else:
        raise ValueError(f"Can't convert to protobuf: {value}")
