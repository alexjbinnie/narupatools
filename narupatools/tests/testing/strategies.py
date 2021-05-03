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

import sys

import numpy as np
from hypothesis import strategies as st

MAX_DOUBLE = sys.float_info.max
MIN_DOUBLE = sys.float_info.min
MAX_FLOAT32 = np.finfo(np.float32).max
MIN_FLOAT32 = np.finfo(np.float32).min
MAX_UINT32 = np.iinfo(np.uint32).max
MAX_INT32 = np.iinfo(np.int32).max
MIN_INT32 = np.iinfo(np.int32).min


def keys():
    return st.text(min_size=1, max_size=10)


def serializable_primitives():
    return st.one_of(
        st.text(),
        st.booleans(),
        st.none(),
        st.floats(min_value=float(MIN_DOUBLE), max_value=float(MAX_DOUBLE)),
        st.integers(min_value=int(MIN_INT32), max_value=int(MAX_INT32)),
    )


def serializable_dictionaries():
    return st.dictionaries(keys(), serializable(), max_size=3)


def serializable_lists():
    return st.lists(serializable(), max_size=3)


def serializable_flat():
    return st.one_of(
        serializable_primitives(), serializable_dictionaries(), serializable_lists()
    )


def serializable():
    return st.recursive(
        serializable_primitives(),
        lambda s: st.one_of(
            st.lists(s, max_size=3), st.dictionaries(keys(), s, max_size=3)
        ),
        max_leaves=3,
    )
