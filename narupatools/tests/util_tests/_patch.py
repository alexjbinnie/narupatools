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

import narupatools.util.patch as patch


def test_patch():
    class Object:
        def func(self) -> str:
            return "A"

    unmodified = Object()
    composed = Object()
    replaced = Object()

    @patch.compose(composed.func)
    def my_func2(self: Object, arg: str) -> str:
        return arg + "B"

    @patch.replace(replaced.func)
    def my_func3(self: Object) -> str:
        return "B"

    assert unmodified.func() == "A"
    assert composed.func() == "AB"
    assert replaced.func() == "B"
