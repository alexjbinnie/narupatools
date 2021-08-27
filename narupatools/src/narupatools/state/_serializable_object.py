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
Defines an object that can be serialized to and from protobuf.

Defines objects that can be easily serialized and deserialized to a form that is
compatible with serialization formats such as JSON or protobuf.

Most serialization methods only work with primitive types (strings, floats and
booleans), along with nest lists and string-valued dictionaries of these types. It is
therefore a common problem of wanting to have a class be able to be serialized into
these formats. This often means writing custom serialization and deserialization code
for each class you define.

SerializableObject is a metaclass for any object which is saying to Narupa that it can
be converted to and from a form that is serializable, and provides the load() and save()
methods to do so.
"""
from abc import abstractmethod
from typing import Protocol, Type, TypeVar, runtime_checkable

from narupatools.state.typing import Serializable

TClass = TypeVar("TClass")


@runtime_checkable
class SerializableObject(Protocol):
    """
    Protocol for objects which can be converted to and from a serializable form.

    Types that are serializable by protobuf consist of primitive types, lists and
    dictionaries.
    """

    @classmethod
    @abstractmethod
    def deserialize(cls: Type[TClass], value: Serializable) -> TClass:
        """
        Deserialize a given object to an instance of this class.

        Keys which share a name with a property will use that property's setter, whilst
        other keys will be stored as arbitrary data which can be accessed using an item
        accessor.

        :param value: Object to deserialize to the given object.
        """
        raise NotImplementedError()

    @abstractmethod
    def serialize(self) -> Serializable:
        """
        Serialize this object to a form that is serializable to protobuf.

        All properties and arbitrary data will be writen to the dictionary.
        """
        raise NotImplementedError()
