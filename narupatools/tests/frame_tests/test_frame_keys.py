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

from narupatools.frame.fields import BondPairs


def test_str_equality():
    assert BondPairs == "bond.pairs"


def test_str_inequality():
    assert BondPairs != "particle.positions"


def test_str_set_membership():
    assert BondPairs in {"bond.pairs"}


def test_key_set_membership():
    assert BondPairs in {BondPairs}


def test_key_str_set_membership():
    assert "bond.pairs" in {BondPairs}


def test_key_dict_key():
    dictionary = {"bond.pairs": 2}
    assert dictionary[BondPairs] == 2
