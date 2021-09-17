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
from typing import Generator

import pytest
from fixtures import SHARED_STATE_WRAPPERS

from narupatools.state import (
    SharedStateCollectionView,
    SharedStateDictionaryView,
    SharedStateReference,
)
from narupatools.state.typing import SerializableDictionary

KEY_OF_STRING_VALUE = "item.abc"
KEY_OF_STRING_VALUE_SHORTHAND = "abc"
INITIAL_STRING_VALUE = "value"
NEW_STRING_VALUE = "new"

KEY_WITH_SAME_STRING_VALUE = "item.alt"

NEW_KEY = "item.def"

KEY_OF_DICT_VALUE = "item.ghi"
KEY_OF_DICT_VALUE_SHORTHAND = "ghi"

DICT_KEY = "item"
DICT_NEW_KEY = "item2"
DICT_INITIAL_VALUE = "value"
DICT_NEW_VALUE = "new"

DICT_VALUE = {DICT_KEY: DICT_INITIAL_VALUE}

NEW_DICT_VALUE = {DICT_KEY: DICT_NEW_VALUE}

ABSENT_KEY = "item.xyz"

KEY_NOT_IN_COLLECTION = "not.an_item"

FULL_DICTIONARY = {
    KEY_OF_STRING_VALUE: INITIAL_STRING_VALUE,
    KEY_WITH_SAME_STRING_VALUE: INITIAL_STRING_VALUE,
    KEY_OF_DICT_VALUE: DICT_VALUE,
    KEY_NOT_IN_COLLECTION: INITIAL_STRING_VALUE,
}

COLLECTION_DICTIONARY = {
    KEY_OF_STRING_VALUE: INITIAL_STRING_VALUE,
    KEY_WITH_SAME_STRING_VALUE: INITIAL_STRING_VALUE,
    KEY_OF_DICT_VALUE: DICT_VALUE,
}


@pytest.fixture(params=SHARED_STATE_WRAPPERS)
def serializable_dictionary(request) -> Generator[SerializableDictionary, None, None]:
    with request.param(copy.deepcopy(FULL_DICTIONARY)) as shared_state:
        yield shared_state


def _create_dict_view(serializable_dictionary):
    return SharedStateDictionaryView(serializable_dictionary)


def _create_collection_view(serializable_dictionary):
    return SharedStateCollectionView(serializable_dictionary, "item.")


@pytest.fixture(params=[_create_dict_view, _create_collection_view])
def view(request, serializable_dictionary):
    return request.param(serializable_dictionary)


@pytest.fixture
def view_dictionary(serializable_dictionary):
    return _create_dict_view(serializable_dictionary)


@pytest.fixture
def view_collection(serializable_dictionary):
    return _create_collection_view(serializable_dictionary)


# SharedStateView.__getitem__ returns the correct SharedStateReference
def test_getitem(serializable_dictionary, view):
    assert view[KEY_OF_STRING_VALUE] == SharedStateReference(
        serializable_dictionary, KEY_OF_STRING_VALUE
    )


# SharedStateView.__getitem__ returns the correct SharedStateReference
def test_getitem_collection_shorthand(serializable_dictionary, view_collection):
    assert view_collection[KEY_OF_STRING_VALUE_SHORTHAND] == SharedStateReference(
        serializable_dictionary, KEY_OF_STRING_VALUE
    )


# SharedStateView.__getitem__ returns a reference where has_value() is True
def test_getitem_reference_hasvalue(serializable_dictionary, view):
    assert view[KEY_OF_STRING_VALUE].has_value()


# SharedStateView.__getitem__ returns a SharedStateReference which in turn has the
# correct snapshot
def test_view_getitem_snapshot(serializable_dictionary, view):
    assert (
        view[KEY_OF_STRING_VALUE].snapshot()
        == serializable_dictionary[KEY_OF_STRING_VALUE]
    )


# SharedStateView.__getitem__ returns a SharedStateReference when a key is absent, and
# the reference does not have a value (has_value() is false)
def test_view_getitem_missingkey_noerror(serializable_dictionary, view):
    reference = view[ABSENT_KEY]
    assert not reference.has_value()


# SharedStateView.__getitem__ throws a KeyError when the key is None
def test_getitem_none_keyerror(serializable_dictionary, view):
    with pytest.raises(KeyError):
        assert view[None]


# SharedStateView.__getitem__ throws a KeyError when the key is not a string
def test_getitem_wrongkeytype_keyerror(serializable_dictionary, view):
    with pytest.raises(KeyError):
        assert view[4]


# SharedStateView.__setitem__ sets a value correctly
def test_setitem(serializable_dictionary, view):
    view[KEY_OF_STRING_VALUE] = NEW_STRING_VALUE
    assert serializable_dictionary[KEY_OF_STRING_VALUE] == NEW_STRING_VALUE


# SharedStateCollectionView.__setitem__ sets a value correctly using a shorthand key
def test_setitem_collection_shorthand(serializable_dictionary, view_collection):
    view_collection[KEY_OF_STRING_VALUE_SHORTHAND] = NEW_STRING_VALUE
    assert serializable_dictionary[KEY_OF_STRING_VALUE] == NEW_STRING_VALUE


# SharedStateView.__setitem__ to update a key affects a reference to that key
def test_setitem_reference(serializable_dictionary, view):
    reference = view[KEY_OF_STRING_VALUE]
    view[KEY_OF_STRING_VALUE] = NEW_STRING_VALUE
    assert reference.has_value()
    assert reference.snapshot() == NEW_STRING_VALUE


# SharedStateView.__setitem__ creates a value if it was not present
def test_setitem_absentkey(serializable_dictionary, view):
    view[NEW_KEY] = NEW_STRING_VALUE
    assert serializable_dictionary[NEW_KEY] == NEW_STRING_VALUE


# SharedStateView.__setitem__ for an absent key affects a reference to that key
def test_setitem_absentkey_reference(serializable_dictionary, view):
    reference = view[NEW_KEY]
    assert not reference.has_value()
    view[NEW_KEY] = NEW_STRING_VALUE
    assert reference.has_value()
    assert reference.snapshot() == NEW_STRING_VALUE


# SharedStateView.__setitem__ cannot be used with a None key
def test_setitem_nonekey(serializable_dictionary, view):
    with pytest.raises(KeyError):
        view[None] = NEW_STRING_VALUE


# SharedStateView.__setitem__ cannot be used with a non-string key
def test_setitem_wrongkeytype(serializable_dictionary, view):
    with pytest.raises(KeyError):
        view[23] = NEW_STRING_VALUE


# SharedStateView.__delitem__ removes value from shared dictionary
def test_delitem(serializable_dictionary, view):
    del view[KEY_OF_STRING_VALUE]
    assert KEY_OF_STRING_VALUE not in view
    assert KEY_OF_STRING_VALUE not in serializable_dictionary


# SharedStateCollectionView.__delitem__ resolves key correctly
def test_delitem_collection_shorthand(serializable_dictionary, view_collection):
    del view_collection[KEY_OF_STRING_VALUE_SHORTHAND]
    assert KEY_OF_STRING_VALUE not in view_collection
    assert KEY_OF_STRING_VALUE not in serializable_dictionary


# SharedStateView.__delitem__ affects reference to that key
def test_delitem_reference(serializable_dictionary, view):
    reference = view[KEY_OF_STRING_VALUE]
    del view[KEY_OF_STRING_VALUE]
    assert not reference.has_value()


# SharedStateView.__delitem__ throws a KeyError if the key is absent
def test_delitem_absentkey(serializable_dictionary, view):
    with pytest.raises(KeyError):
        del view[ABSENT_KEY]


# SharedStateView.__delitem__ throws a KeyError if the key is None
def test_delitem_nonekey(serializable_dictionary, view):
    with pytest.raises(KeyError):
        del view[None]


# SharedStateView.__delitem__ throws a KeyError if the key is not a string
def test_delitem_wrongkeytype(serializable_dictionary, view):
    with pytest.raises(KeyError):
        del view[42]


# SharedStateView.set() overrides an existing item
def test_set(serializable_dictionary, view):
    view.set(KEY_OF_STRING_VALUE, NEW_STRING_VALUE)
    assert serializable_dictionary[KEY_OF_STRING_VALUE] == NEW_STRING_VALUE


# SharedStateCollectionView.set() resolves the key correctly
def test_set_collection_shorthand(serializable_dictionary, view_collection):
    view_collection.set(KEY_OF_STRING_VALUE_SHORTHAND, NEW_STRING_VALUE)
    assert serializable_dictionary[KEY_OF_STRING_VALUE] == NEW_STRING_VALUE


# SharedStateView.set() returns a correct reference to the item
def test_set_reference(serializable_dictionary, view):
    reference = view.set(KEY_OF_STRING_VALUE, NEW_STRING_VALUE)
    assert reference.key == KEY_OF_STRING_VALUE
    assert reference.has_value()
    assert reference.snapshot() == NEW_STRING_VALUE


# SharedStateView.set() returns a reference that equals one obtained before set was
# called
def test_set_reference_equality(serializable_dictionary, view):
    old_reference = view[KEY_OF_STRING_VALUE]
    reference = view.set(KEY_OF_STRING_VALUE, NEW_STRING_VALUE)
    assert old_reference == reference


# SharedStateView.set() affects an existing reference
def test_set_absentkey_reference(serializable_dictionary, view):
    reference = view[ABSENT_KEY]
    view.set(ABSENT_KEY, NEW_STRING_VALUE)
    assert reference.has_value()
    assert reference.snapshot() == NEW_STRING_VALUE


# SharedStateView.update() modifies a value in an existing dict
def test_update_dict(serializable_dictionary, view):
    view.update(KEY_OF_DICT_VALUE, **{DICT_KEY: DICT_NEW_VALUE})
    assert serializable_dictionary[KEY_OF_DICT_VALUE] == {DICT_KEY: DICT_NEW_VALUE}


# SharedStateCollectionView.update() resolves keys correctly
def test_update_collection_shorthand(serializable_dictionary, view_collection):
    view_collection.update(KEY_OF_DICT_VALUE_SHORTHAND, **{DICT_KEY: DICT_NEW_VALUE})
    assert serializable_dictionary[KEY_OF_DICT_VALUE] == {DICT_KEY: DICT_NEW_VALUE}


# SharedStateView.update() adds a new value to an existing dict
def test_update_dict_newkey(serializable_dictionary, view):
    view.update(KEY_OF_DICT_VALUE, **{DICT_NEW_KEY: DICT_NEW_VALUE})
    assert serializable_dictionary[KEY_OF_DICT_VALUE] == {
        DICT_KEY: DICT_INITIAL_VALUE,
        DICT_NEW_KEY: DICT_NEW_VALUE,
    }


# SharedStateView.update() throws a TypeError if the value does not support item
# assignment
def test_update_cantdoitemassignment(serializable_dictionary, view):
    with pytest.raises(TypeError):
        view.update(KEY_OF_STRING_VALUE, **{DICT_NEW_KEY: DICT_NEW_VALUE})


# SharedStateView.update() throws a TypeError if the value does not support item
# assignment
def test_update_notmutablemapping(serializable_dictionary, view):
    with pytest.raises(TypeError):
        view.update(KEY_OF_STRING_VALUE, **{DICT_NEW_KEY: DICT_NEW_VALUE})


# SharedStateView.update() throws a KeyError if the key is not present
def test_update_absentkey(serializable_dictionary, view):
    with pytest.raises(KeyError):
        view.update(ABSENT_KEY, **{DICT_NEW_KEY: DICT_NEW_VALUE})


# SharedStateView.update() throws a KeyError if the key is None
def test_update_nonekey(serializable_dictionary, view):
    with pytest.raises(KeyError):
        view.update(None, **{DICT_NEW_KEY: DICT_NEW_VALUE})


# SharedStateView.update() throws a KeyError if the key is not a string
def test_update_wrongkeytype(serializable_dictionary, view):
    with pytest.raises(KeyError):
        view.update(223, **{DICT_NEW_KEY: DICT_NEW_VALUE})


# SharedStateView.update() works with no keyword arguments
def test_update_nokwargs(serializable_dictionary, view):
    view.update(KEY_OF_DICT_VALUE)
    assert serializable_dictionary[KEY_OF_DICT_VALUE] == DICT_VALUE


# SharedStateView.modify() returns a snapshot of the value
def test_modify_string(serializable_dictionary, view):
    with view.modify(KEY_OF_STRING_VALUE) as value:
        assert value == INITIAL_STRING_VALUE


# SharedStateCollectionView.modify() resolves keys correctly
def test_modify_collection_shorthand(serializable_dictionary, view_collection):
    with view_collection.modify(KEY_OF_STRING_VALUE_SHORTHAND) as value:
        assert value == INITIAL_STRING_VALUE


# SharedStateView.modify() applies modifications to dictionaries
def test_modify_dict(serializable_dictionary, view):
    with view.modify(KEY_OF_DICT_VALUE) as value:
        assert value == {DICT_KEY: DICT_INITIAL_VALUE}
        value[DICT_KEY] = DICT_NEW_VALUE
    assert serializable_dictionary[KEY_OF_DICT_VALUE] == {DICT_KEY: DICT_NEW_VALUE}


# SharedStateView.modify() affects a previously defined reference
def test_modify_dict_reference(serializable_dictionary, view):
    reference = view[KEY_OF_DICT_VALUE]
    with view.modify(KEY_OF_DICT_VALUE) as value:
        assert value == {DICT_KEY: DICT_INITIAL_VALUE}
        value[DICT_KEY] = DICT_NEW_VALUE
    assert reference.snapshot() == {DICT_KEY: DICT_NEW_VALUE}


# SharedStateView.modify() throws a KeyError for a missing key
def test_modify_absentkey(serializable_dictionary, view):
    with pytest.raises(KeyError), view.modify(ABSENT_KEY):
        pass


# SharedStateView.modify() throws a KeyError for a key of None
def test_modify_nonekey(serializable_dictionary, view):
    with pytest.raises(KeyError), view.modify(None):
        pass


# SharedStateView.modify() throws a KeyError for a key which is not a string
def test_modify_wrongkeytype(serializable_dictionary, view):
    with pytest.raises(KeyError), view.modify(231):
        pass


# SharedStateView.snapshot() returns a snapshot of the whole dictionary
def test_snapshot(serializable_dictionary, view_dictionary):
    snapshot = view_dictionary.snapshot()
    assert serializable_dictionary == snapshot


# SharedStateCollectionView.snapshot() returns a snapshot of just the collection
def test_snapshot_collection(serializable_dictionary, view_collection):
    snapshot = view_collection.snapshot()
    assert COLLECTION_DICTIONARY == snapshot


# SharedStateReference.has_value() is true if the key is present
def test_reference_has_value(serializable_dictionary, view):
    reference = view[KEY_OF_STRING_VALUE]
    assert reference.has_value()


# SharedStateReference.has_value() is false if the key is not present
def test_reference_has_no_value(serializable_dictionary, view):
    reference = view[ABSENT_KEY]
    assert not reference.has_value()


# SharedStateReference.snapshot() returns current value
def test_reference_snapshot(serializable_dictionary, view):
    reference = view[KEY_OF_STRING_VALUE]
    assert reference.snapshot() == INITIAL_STRING_VALUE


# SharedStateReference.snapshot() returns a copy of its value, so modifications of the
# dictionary don't affect its value
def test_reference_snapshot_iscopy(serializable_dictionary, view):
    reference = view[KEY_OF_STRING_VALUE]
    value = reference.snapshot()
    reference.set(NEW_STRING_VALUE)
    assert value == INITIAL_STRING_VALUE


# SharedStateReference.snapshot() returns a copy of its value, so modifications to it
# don't affect the dictionary
def test_reference_snapshot_iscopy_dict(serializable_dictionary, view):
    reference = view[KEY_OF_DICT_VALUE]
    value = reference.snapshot()
    value[DICT_KEY] = DICT_NEW_VALUE
    assert serializable_dictionary[KEY_OF_DICT_VALUE] == {DICT_KEY: DICT_INITIAL_VALUE}


# SharedStateReference.remove() removes the value
def test_reference_remove(serializable_dictionary, view):
    reference = view[KEY_OF_STRING_VALUE]
    reference.remove()
    assert KEY_OF_STRING_VALUE not in serializable_dictionary


# SharedStateReference.remove() raises a KeyError if there isn't a value
def test_reference_remove_absent(serializable_dictionary, view):
    reference = view[ABSENT_KEY]
    with pytest.raises(KeyError):
        reference.remove()


# SharedStateReference.set() overrides an new value
def test_reference_set(serializable_dictionary, view):
    reference = view[KEY_OF_STRING_VALUE]
    reference.set(NEW_STRING_VALUE)
    assert serializable_dictionary[KEY_OF_STRING_VALUE] == NEW_STRING_VALUE


# SharedStateReference.set() sets a new value even when there isn't a value already
def test_reference_set_absent(serializable_dictionary, view):
    reference = view[ABSENT_KEY]
    reference.set(NEW_STRING_VALUE)
    assert serializable_dictionary[ABSENT_KEY] == NEW_STRING_VALUE


# SharedStateReference.modify() returns a dict that can be edited
def test_reference_modify_dict(serializable_dictionary, view):
    reference = view[KEY_OF_DICT_VALUE]
    with reference.modify() as dict_:
        dict_[DICT_KEY] = DICT_NEW_VALUE
    assert serializable_dictionary[KEY_OF_DICT_VALUE] == NEW_DICT_VALUE


# SharedStateReference.modify() throws a KeyError if the reference doesn't have a value
def test_reference_modify_absentkey(serializable_dictionary, view):
    reference = view[ABSENT_KEY]
    with pytest.raises(KeyError), reference.modify():
        pass


# SharedStateReference.update() can modify existing keys
def test_reference_update_dict_replacekey(serializable_dictionary, view):
    reference = view[KEY_OF_DICT_VALUE]
    reference.update(**{DICT_KEY: DICT_NEW_VALUE})
    assert serializable_dictionary[KEY_OF_DICT_VALUE] == {DICT_KEY: DICT_NEW_VALUE}


# SharedStateReference.update() can add new keys to a dictionary value
def test_reference_update_dict_newkey(serializable_dictionary, view):
    reference = view[KEY_OF_DICT_VALUE]
    reference.update(**{DICT_NEW_KEY: DICT_NEW_VALUE})
    assert serializable_dictionary[KEY_OF_DICT_VALUE] == {
        DICT_KEY: DICT_INITIAL_VALUE,
        DICT_NEW_KEY: DICT_NEW_VALUE,
    }


# SharedStateReference.update() throws a TypeError if the value does not support item
# assignment
def test_reference_update_string(serializable_dictionary, view):
    reference = view[KEY_OF_STRING_VALUE]
    with pytest.raises(TypeError):
        reference.update(key="new")


# SharedStateReference.update() throws a KeyError if there is no value
def test_reference_update_absentkey(serializable_dictionary, view):
    reference = view[ABSENT_KEY]
    with pytest.raises(KeyError):
        reference.update(key="new")


# SharedStateReference supports equality
def test_reference_equality(serializable_dictionary, view):
    assert view[KEY_OF_STRING_VALUE] == view[KEY_OF_STRING_VALUE]


# SharedStateReference supports inequality, where references have the same value but
# different keys
def test_reference_inequality(serializable_dictionary, view):
    assert view[KEY_OF_STRING_VALUE] != view[KEY_WITH_SAME_STRING_VALUE]
    assert (
        serializable_dictionary[KEY_OF_STRING_VALUE]
        == serializable_dictionary[KEY_WITH_SAME_STRING_VALUE]
    )


# SharedStateView.__contains__ returns true if the given reference's key is in the
# dictionary
def test_reference_in_collection(serializable_dictionary, view):
    assert view[KEY_OF_STRING_VALUE] in view


# SharedStateView.__contains__ returns false if the given reference is no longer present
def test_reference_not_in_collection(serializable_dictionary, view):
    item = view[KEY_OF_STRING_VALUE]
    del view[KEY_OF_STRING_VALUE]
    assert item not in view


# SharedStateReference.__repr__ returns the expected result
def test_reference_repr(serializable_dictionary, view):
    reference = view[KEY_OF_STRING_VALUE]
    assert repr(reference) == f"<SharedStateReference key:{KEY_OF_STRING_VALUE}>"


# SharedStateView.__contains__ returns true is the string key is present
def test_key_in_collection(serializable_dictionary, view):
    assert KEY_OF_STRING_VALUE in view


# SharedStateCollectionView.__contains__ returns true is the shorthand key is present
def test_key_in_collection_shorthand(serializable_dictionary, view_collection):
    assert KEY_OF_STRING_VALUE_SHORTHAND in view_collection


# SharedStateView.__contains__ returns false if the string key is not present
def test_key_not_in_collection(serializable_dictionary, view):
    assert ABSENT_KEY not in view


# SharedStateView.__contains__ returns false if the string key was deleted
def test_deleted_key_not_in_collection(serializable_dictionary, view):
    del view[KEY_OF_STRING_VALUE]
    assert KEY_OF_STRING_VALUE not in view


# SharedStateView.__contains__ returns false for a key of None
def test_nonekey_not_in_collection(serializable_dictionary, view):
    assert None not in view


# SharedStateView.__contains__ returns false for a key of a non-string type
def test_wrongkeytype_not_in_collection(serializable_dictionary, view):
    assert 231 not in view


# SharedStateView.__len__ returns the correct length
def test_len(serializable_dictionary, view_dictionary):
    assert len(FULL_DICTIONARY) == len(view_dictionary)


# SharedStateCollectionView.__len__ returns the number of keys with the correct prefix
def test_len_collection(serializable_dictionary, view_collection):
    assert len(COLLECTION_DICTIONARY) == len(view_collection)


# SharedStateView.clear removes all values
def test_clear(serializable_dictionary, view_dictionary):
    view_dictionary.clear()
    assert len(serializable_dictionary) == 0
    assert len(view_dictionary) == 0


# SharedStateView.clear removes all values
def test_clear_collection(serializable_dictionary, view_collection):
    view_collection.clear()
    assert len(serializable_dictionary) == len(FULL_DICTIONARY) - len(
        COLLECTION_DICTIONARY
    )
    assert len(view_collection) == 0


# SharedStateDictionaryView implements keys() correctly
def test_keys(serializable_dictionary, view_dictionary):
    assert view_dictionary.keys() == FULL_DICTIONARY.keys()


# SharedStateCollectionView implements keys() correctly
def test_keys_collection(serializable_dictionary, view_collection):
    assert view_collection.keys() == COLLECTION_DICTIONARY.keys()


# SharedStateDictionaryView implements values() correctly
def test_values(serializable_dictionary, view_dictionary):
    assert set(view_dictionary.values()) == {
        view_dictionary[key] for key in FULL_DICTIONARY.keys()
    }


# SharedStateCollectionView implements values() correctly
def test_values_collection(serializable_dictionary, view_collection):
    assert set(view_collection.values()) == {
        view_collection[key] for key in COLLECTION_DICTIONARY.keys()
    }


# SharedStateDictionaryView implements items() correctly
def test_items(serializable_dictionary, view_dictionary):
    assert dict(view_dictionary.items()) == {
        key: view_dictionary[key] for key in FULL_DICTIONARY.keys()
    }


# SharedStateCollectionView implements items() correctly
def test_items_collection(serializable_dictionary, view_collection):
    assert dict(view_collection.items()) == {
        key: view_collection[key] for key in COLLECTION_DICTIONARY.keys()
    }


# SharedStateView implements get() correctly
def test_get(serializable_dictionary, view):
    reference = view.get(KEY_OF_STRING_VALUE)
    assert reference.snapshot() == INITIAL_STRING_VALUE


@pytest.mark.parametrize(
    ("test_input", "expected"),
    [
        ({}, "<SharedStateDictionaryView 0 item(s)>"),
        ({"abc": 2}, "<SharedStateDictionaryView 1 item(s)>"),
        ({"abc": 1, "def": 3}, "<SharedStateDictionaryView 2 item(s)>"),
        (
            {"abc": True, "def": -1.2, "ghi": None},
            "<SharedStateDictionaryView 3 item(s)>",
        ),
    ],
)
def test_repr(test_input, expected):
    assert repr(SharedStateDictionaryView(test_input)) == expected


@pytest.mark.parametrize(
    ("test_input", "expected"),
    [
        ({}, "<SharedStateCollectionView 0 item(s)>"),
        ({"item.abc": 2}, "<SharedStateCollectionView 1 item(s)>"),
        ({"item.abc": 1, "item.def": 3}, "<SharedStateCollectionView 2 item(s)>"),
        (
            {"abc": True, "item.def": -1.2, "ghi": None},
            "<SharedStateCollectionView 1 item(s)>",
        ),
    ],
)
def test_repr_collection(test_input, expected):
    assert repr(SharedStateCollectionView(test_input, prefix="item.")) == expected
