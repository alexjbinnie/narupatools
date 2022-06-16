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

from narupa.trajectory import FrameData

from narupatools.frame import frame_to_pdb_string


def test_pdb_empty_frame():
    frame_to_pdb_string(FrameData())


def test_pdb_line_length(frame):
    pdb = frame_to_pdb_string(frame)
    lines = pdb.splitlines()
    for line in lines:
        assert len(line) == 80, f"Line is not of length 80: `{line}`"
