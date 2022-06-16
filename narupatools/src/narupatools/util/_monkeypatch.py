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

from typing import Callable, Type


def monkeypatch(baseclass: Type) -> Callable[[Type], None]:
    """Decorate that modifies a base class via monkeypatching."""

    def patch(derivedclass: Type) -> None:
        for clss in derivedclass.mro():
            if clss == baseclass or clss in baseclass.mro():
                continue
            for key, value in clss.__dict__.items():
                if key in {"__doc__", "__module__", "__dict__", "__weakref__"}:
                    continue
                setattr(baseclass, key, value)

    return patch
