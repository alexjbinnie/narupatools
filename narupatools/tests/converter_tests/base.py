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

from abc import ABCMeta, abstractmethod
from typing import AbstractSet

import numpy as np
import pytest

from narupatools.frame import FrameData, convert
from narupatools.frame.fields import BondPairs, FrameKey


def check_equal(field, original, new):
    if property is BondPairs:
        assert np.array_equiv(np.sort(original, axis=0), np.sort(new, axis=0))
    elif isinstance(original, float) or (
        isinstance(original, np.ndarray) and original.dtype == float
    ):
        assert new == pytest.approx(original)
    elif isinstance(original, np.ndarray):
        assert np.all(new == original)
    else:
        assert new == original


class BaseTestConverter(metaclass=ABCMeta):
    frame_to_object_fields: AbstractSet[FrameKey] = set()

    object_to_frame_fields: AbstractSet[FrameKey] = set()

    minimal_fields: AbstractSet[FrameKey] = set()

    @property
    @abstractmethod
    def object_type(self):
        pass

    def test_round_trip(self, frame):
        topology = convert(frame, self.object_type)
        frame_conv = convert(topology, FrameData)
        for field in self.frame_to_object_fields & self.object_to_frame_fields:
            assert field.key in frame_conv

            original = field.get(frame)
            new = field.get(frame_conv)

            check_equal(field, original, new)

    def test_minimal(self, make_frame):
        frame = make_frame(self.minimal_fields)
        topology = convert(frame, self.object_type)
        frame_conv = convert(topology, FrameData)
        for field in self.minimal_fields:
            assert field.key in frame_conv

            original = field.get(frame)
            new = field.get(frame_conv)

            check_equal(field, original, new)
