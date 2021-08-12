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

_CATEGORIES = [
    "atom",
    "integrate",
    "minimize",
    "pair",
    "bond",
    "angle",
    "dihedral",
    "improper",
    "kspace",
    "compute",
    "fix",
    "region",
    "dump",
    "command",
]

_VALID_OPTIONS = [
    "full",
    "verlet",
    "cg",
    "lj/cut",
    "harmonic",
    "charmm",
    "class2",
    "cossq",
    "ewald",
    "angle",
    "external",
    "cone",
    "atom",
    "read_data",
]


@pytest.mark.parametrize("category", _CATEGORIES)
def test_has_styles(lammps, category):
    assert len(getattr(lammps, f"{category}_styles")) > 0


@pytest.mark.parametrize(
    ("category", "valid_style"), list(zip(_CATEGORIES, _VALID_OPTIONS))
)
def test_valid_style(lammps, category, valid_style):
    assert valid_style in getattr(lammps, f"{category}_styles")


@pytest.mark.parametrize("category", _CATEGORIES)
def test_invalid_style(lammps, category):
    assert "invalid_style" not in getattr(lammps, f"{category}_styles")


@pytest.mark.parametrize("category", _CATEGORIES)
def test_style_index(lammps, category):
    assert isinstance(getattr(lammps, f"{category}_styles")[0], str)


@pytest.mark.parametrize("category", _CATEGORIES)
def test_style_negative_index(lammps, category):
    with pytest.raises(IndexError):
        _ = getattr(lammps, f"{category}_styles")[-1]


@pytest.mark.parametrize("category", _CATEGORIES)
def test_style_outofrange_index(lammps, category):
    with pytest.raises(IndexError):
        _ = getattr(lammps, f"{category}_styles")[
            len(getattr(lammps, f"{category}_styles"))
        ]
