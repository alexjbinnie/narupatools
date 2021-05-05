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

from typing import Iterable, Tuple

import numpy as np

class Element:
    number: int

class Atom:
    index: int
    residue: Residue
    element: Element
    name: str

class Residue:
    name: str
    chain: Chain
    index: int

class Chain:
    index: int

class Topology:
    chains: Iterable[Chain]
    residues: Iterable[Residue]
    atoms: Iterable[Atom]
    bonds: Iterable[Tuple[Atom, Atom]]
    @property
    def n_atoms(self) -> int: ...
    @property
    def n_residues(self) -> int: ...
    @property
    def n_chains(self) -> int: ...
    @property
    def n_bonds(self) -> int: ...

class Trajectory:
    xyz: np.ndarray
    n_atoms: int
    n_residues: int
    n_chains: int
    n_bonds: int
    topology: Topology
    timestep: float
    unitcell_vectors: np.ndarray
    def __len__(self) -> int: ...

def load(file: str) -> Trajectory: ...
