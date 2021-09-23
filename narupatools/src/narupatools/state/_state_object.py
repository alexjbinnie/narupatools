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
Generic implementation of a serializable object.

SharedStateObject is an implementation of a SerializableObject, that can store arbitrary
data as a dictionary, but also automatically use any properties defined in the class. In
most cases, this class is the appropriate one to subclass for defining a specific kind
of object that's storable in the shared state dictionary.

.. code:: python

    from narupa.state.state_object import SharedStateObject


    # Define a custom class, with a property x that must be non-negative.
    class MyObject(SharedStateObject):

        _x: int

        def __init__(self, **kwargs):
            super().__init__(**kwargs)
            x = 0

        @property
        def x(self) -> int:
            return self._x

        @x.setter
        def x(self, value):
            if value < 0:
                value = 0
            self._x = value


    # Define a dictionary we want to deserialize
    dictionary = {'x': -1, 'y': 2}

    # Create the object using the class method load()
    myobj = MyObject.load(dictionary)

    # Check that the property setter was called
    print(myobj['x'])  # we could also use myobj.x

    # Check that y was also stored
    print(myobj['y'])
"""
import contextlib
from collections.abc import Mapping
from typing import Any, ClassVar, Dict, Union

from ..override import override
from .typing import Serializable, SerializableObject


class SharedStateObject(SerializableObject):
    """
    Base class for an object which is stored in the shared state dictionary.

    This stores arbitrary fields of serializable data. Any properties that are defined
    on a subclass of SharedStateObject are serialized automatically, calling their get
    and set methods as appropriate.
    """

    _arbitrary_data: Dict[str, Serializable]
    _serializable_properties: ClassVar[Dict[str, property]] = {}

    def __init__(self, **kwargs: Serializable):
        self._arbitrary_data = {}
        for key, value in kwargs.items():
            if key in self._serializable_properties:
                self._serializable_properties[key].fset(self, value)  # type: ignore[misc]
            else:
                self._arbitrary_data[key] = value

    def __init_subclass__(cls, **kwargs: Any):
        # Find all properties in this class
        cls._serializable_properties = {}
        for key in dir(cls):
            if key.startswith("_"):
                continue
            value = getattr(cls, key)
            if isinstance(value, property):
                cls._serializable_properties[key] = value

    @classmethod
    @override
    def deserialize(cls, value: Serializable) -> "SharedStateObject":  # noqa: D102
        if isinstance(value, Mapping):
            return cls(**value)
        raise ValueError

    @override
    def serialize(self) -> Dict[str, Serializable]:  # noqa: D102
        dictionary = {k: v for k, v in self._arbitrary_data.items()}
        for key, python_property in self.__class__._serializable_properties.items():
            with contextlib.suppress(AttributeError):
                value = python_property.fget(self)  # type: ignore
                if value is not None:  # Only save non None values
                    dictionary[key] = value

        return dictionary

    def __setitem__(self, key: str, value: Union[Serializable, Any]) -> None:
        """
        Set a value for this object.

        This uses a property if the key corresponds to a property of this object and
        using an internal dictionary for other arbitrary values

        :param key: The key to store the value under. This will be the key used in the
                    serialized form of this object
        :param value: The value to be stored. If key does not correspond to a property,
                      this value must be serializable. If the key does correspond to a
                      property, then the properties setter must be able to handle
                      converting the value to something serializable
        """
        if key in self.__class__._serializable_properties:
            self.__class__._serializable_properties[key].fset(
                self, value
            )  # type: ignore
        else:
            self._arbitrary_data[key] = value

    def __getitem__(self, key: str) -> Serializable:
        """
        Get a value for this object.

        This uses a property if the key corresponds to a property of this object and
        using an internal dictionary for other arbitrary values

        :param key: The key to store the value under. This will be the key used in the
                    serialized form of this object

        :raises KeyError: If the given key is not present
        """
        if key in self.__class__._serializable_properties:
            return self.__class__._serializable_properties[key].fget(  # type: ignore
                self
            )
        else:
            return self._arbitrary_data[key]

    def __getattr__(self, key: str) -> Serializable:
        try:
            return self._arbitrary_data[key]
        except KeyError as e:
            raise AttributeError from e

    def __eq__(self, other: Any) -> bool:
        return (
            isinstance(other, SharedStateObject)
            and self._arbitrary_data == other._arbitrary_data
        )
