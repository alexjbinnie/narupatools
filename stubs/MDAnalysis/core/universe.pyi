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
from typing import Any, Optional, Sequence, Union

from MDAnalysis import AtomGroup
from MDAnalysis.core.topology import Topology
from MDAnalysis.core.topologyobjects import TopologyGroup

class Universe:
    def __init__(
        self,
        file: Union[str, Topology],
        guess_bonds: bool = ...,
        format: Optional[str] = ...,
    ): ...
    def select_atoms(self, sel: str) -> AtomGroup: ...
    @property
    def dimensions(self) -> Sequence[float]: ...
    @property
    def atoms(self) -> AtomGroup: ...
    @property
    def bonds(self) -> TopologyGroup: ...
    @property
    def trajectory(self) -> Any: ...
    @trajectory.setter
    def trajectory(self, value: Any) -> None: ...
