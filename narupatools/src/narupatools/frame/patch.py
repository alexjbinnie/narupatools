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

"""Patch to FrameData that allows copy() to work with empty arrays."""

from narupa.trajectory import FrameData


def _copy(self: FrameData) -> FrameData:
    frame = FrameData()
    for key, value in self.raw.arrays.items():
        frame.raw.arrays[key].CopyFrom(value)
    for key, value in self.raw.values.items():
        frame.raw.values[key].CopyFrom(value)
    return frame


FrameData.copy = _copy  # type: ignore
