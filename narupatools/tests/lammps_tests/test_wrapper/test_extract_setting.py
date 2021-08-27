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

import pytest

pytest.importorskip("lammps")

from narupatools.lammps.exceptions import SettingNotFoundError


def test_extract_setting(lammps):
    value = lammps.extract_setting("nlocal")
    assert value == 2004


def test_extract_setting_missing(lammps):
    with pytest.raises(SettingNotFoundError):
        lammps.extract_setting("missing_setting")
