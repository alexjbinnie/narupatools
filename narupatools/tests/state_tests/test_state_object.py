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

import copy
from typing import Generator, MutableMapping

import pytest
from fixtures import SHARED_STATE_WRAPPERS
from hypothesis import given
from strategies import keys, serializable, serializable_dictionaries

from narupatools.state import (
    SerializableObject,
    SharedStateDictionaryView,
    SharedStateObject,
)
from narupatools.state.typing import Serializable


class ExampleStateObject(SharedStateObject):
    def __init__(self, named_field=None, **kwargs):
        super().__init__(**kwargs)
        self.named_field = named_field

    @property
    def named_field(self):
        return self._named_field

    @named_field.setter
    def named_field(self, value):
        self._named_field = value

    def __eq__(self, other):
        return (
            isinstance(other, ExampleStateObject)
            and self.named_field == other.named_field
            and self._arbitrary_data == other._arbitrary_data
        )

    def __repr__(self):
        representation = "<ExampleStateObject"
        representation += f" named_field:{self.named_field}"
        if self._arbitrary_data:
            representation += f" extra:{self._arbitrary_data}"
        representation += ">"
        return representation


KEY = "item.item1"
VALUE = {"named_field": "item"}
VALUE_OBJ = ExampleStateObject(named_field="item")

KEY_EXTRA = "item.item2"
VALUE_EXTRA = {"named_field": "item", "extra_data": "extra_data"}
VALUE_EXTRA_OBJ = ExampleStateObject(named_field="item", extra_data="extra_data")

KEY_NON_COLLECTION = "not_in_collection"
VALUE_NON_COLLECTION = {"named_field": "not_in_collection"}
VALUE_NON_COLLECTION_OBJ = ExampleStateObject(named_field="not_in_collection")

KEY_NEW = "item.new"
KEY_NEW_SHORTHAND = "new"
VALUE_NEW = {"named_field": "item new"}
VALUE_NEW_OBJ = ExampleStateObject(named_field="item new")

DICT = {KEY_NON_COLLECTION: VALUE_NON_COLLECTION, KEY: VALUE, KEY_EXTRA: VALUE_EXTRA}


@pytest.fixture
def object() -> ExampleStateObject:
    return ExampleStateObject(named_field="item")


@pytest.fixture(params=SHARED_STATE_WRAPPERS)
def wrapped(request) -> Generator[MutableMapping[str, Serializable], None, None]:
    with request.param(copy.deepcopy(DICT)) as shared_state:
        yield shared_state


@pytest.fixture
def view(wrapped):
    return SharedStateDictionaryView(wrapped)


@pytest.fixture
def collection(view):
    return view.collection("item.", ExampleStateObject)


# Getting a snapshot of a random key returns a dictionary
def test_ref_snapshot_raw(view):
    assert view[KEY_NON_COLLECTION].snapshot() == VALUE_NON_COLLECTION


# Can cast snapshot to arbitrary object using second argument of SharedStateView.get
def test_ref_snapshot_typed(view):
    assert (
        view.get(KEY_NON_COLLECTION, ExampleStateObject).snapshot()
        == VALUE_NON_COLLECTION_OBJ
    )


# Snapshots from references from collections are automatically cast
def test_collection_item_snapshot(collection):
    assert collection[KEY].snapshot() == VALUE_OBJ


# SharedStateCollectionView.snapshot() gives dictionary
def test_collection_full_snapshot(collection):
    assert collection.snapshot() == {KEY: VALUE_OBJ, KEY_EXTRA: VALUE_EXTRA_OBJ}


# SharedStateCollectionView.set() accepts argument that is raw dictionary
def test_collection_set_raw(view, collection):
    collection.set(KEY_NEW, VALUE_NEW)
    assert view[KEY_NEW].snapshot() == VALUE_NEW
    assert collection[KEY_NEW].snapshot() == VALUE_NEW_OBJ


# SharedStateCollectionView.set accepts argument that is correct type
def test_collection_set_typed(view, collection):
    collection.set(KEY_NEW, VALUE_NEW_OBJ)
    assert view[KEY_NEW].snapshot() == VALUE_NEW
    assert collection[KEY_NEW].snapshot() == VALUE_NEW_OBJ


# SharedStateReference.set accepts argument that is raw dictionary
def test_collection_reference_set_raw(view, collection):
    collection[KEY].set(VALUE_NEW)
    assert view[KEY].snapshot() == VALUE_NEW
    assert collection[KEY].snapshot() == VALUE_NEW_OBJ


# SharedStateReference.set accepts argument that is correct type
def test_collection_reference_set_typed(view, collection):
    collection[KEY].set(VALUE_NEW_OBJ)
    assert view[KEY].snapshot() == VALUE_NEW
    assert collection[KEY].snapshot() == VALUE_NEW_OBJ


# SharedStateCollectionView.update() works with SharedStateObject and modifies property
def test_collection_update_property(view, collection):
    collection.update(KEY, **{"named_field": VALUE_NEW["named_field"]})
    assert view[KEY].snapshot() == VALUE_NEW
    assert collection[KEY].snapshot() == VALUE_NEW_OBJ


# SharedStateCollectionView.update() works with SharedStateObject and modifies
# arbitrary data
def test_collection_update_arbitrary(view, collection):
    collection.update(KEY, **{"extra_data": "extra_data"})
    assert view[KEY].snapshot() == VALUE_EXTRA
    assert collection[KEY].snapshot() == VALUE_EXTRA_OBJ


# SharedStateCollectionView.update() works with SharedStateObject and modifies property
def test_collection_reference_update_property(view, collection):
    collection[KEY].update(**{"named_field": VALUE_NEW["named_field"]})
    assert view[KEY].snapshot() == VALUE_NEW
    assert collection[KEY].snapshot() == VALUE_NEW_OBJ


# SharedStateCollectionView.update() works with SharedStateObject and modifies
# arbitrary data
def test_collection_reference_update_arbitrary(view, collection):
    collection[KEY].update(**{"extra_data": "extra_data"})
    assert view[KEY].snapshot() == VALUE_EXTRA
    assert collection[KEY].snapshot() == VALUE_EXTRA_OBJ


# SharedStateCollectionView.add accepts argument that is raw dictionary
def test_collection_add_raw(view, collection):
    collection.add(VALUE_NEW)
    assert VALUE_NEW in view.snapshot().values()
    assert VALUE_NEW_OBJ in collection.snapshot().values()


# SharedStateCollectionView.add accepts argument that is correct type
def test_collection_add_typed(view, collection):
    collection.add(VALUE_NEW_OBJ)
    assert VALUE_NEW in view.snapshot().values()
    assert VALUE_NEW_OBJ in collection.snapshot().values()


# SharedStateReference.modify() works with assigning properties
def test_collection_modify(view, collection):
    with collection[KEY].modify() as item:
        item.named_field = VALUE_NEW["named_field"]
    assert view[KEY].snapshot() == VALUE_NEW
    assert collection[KEY].snapshot() == VALUE_NEW_OBJ


@given(serializable_dictionaries())
def test_load_then_save(dict):
    obj = ExampleStateObject.deserialize(dict)
    assert obj.serialize() == dict


@given(serializable_dictionaries(), keys(), serializable())
def test_arbitrary_data(dict, key, value):
    obj = ExampleStateObject.deserialize(dict)
    dict[key] = value
    obj[key] = value
    new_dict = obj.serialize()
    assert new_dict == dict
    new_obj = ExampleStateObject.deserialize(new_dict)
    assert new_obj == obj


def set_named_field_setitem(object):
    object["named_field"] = "new"
    assert object.named_field == "new"


def get_named_field_getitem(object):
    assert object["named_field"] == object.named_field


def test_ref_repr(collection):
    ref = collection[KEY]
    assert repr(ref) == f"<SharedStateReference key:{KEY} type:ExampleStateObject>"


class UnimplementedSerializable(SerializableObject):
    pass


def test_serializable_unimplemented():
    with pytest.raises(TypeError):
        UnimplementedSerializable()
    with pytest.raises(NotImplementedError):
        UnimplementedSerializable.deserialize({})
