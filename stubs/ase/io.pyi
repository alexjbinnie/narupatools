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

from pathlib import PurePath
from typing import IO, List, Optional, Sequence, Union

from ase.atoms import Atoms

NameOrFile = Union[str, PurePath, IO]

def write(
    filename: NameOrFile, images: Union[Atoms, Sequence[Atoms]], format: str = ...
) -> None: ...
def read(
    filename: NameOrFile, format: Optional[str] = ...
) -> Union[Atoms, List[Atoms]]: ...
