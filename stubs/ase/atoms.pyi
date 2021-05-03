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

from typing import Any, Dict, List, Optional, Sequence, Union

import numpy as np
from ase.calculators.calculator import Calculator

class Atom:
    index: int
    number: int
    position: np.ndarray
    def __init__(self, symbols: str, position: np.ndarray, mass: float): ...

class Atoms:
    constraints: List[Any]
    arrays: Dict[str, Any]
    positions: np.ndarray
    def __init__(
        self,
        symbols: Optional[Union[Sequence[str], np.ndarray]] = ...,
        positions: Optional[np.ndarray] = ...,
        velocities: Optional[np.ndarray] = ...,
        cell: Optional[np.ndarray] = ...,
        charges: Optional[np.ndarray] = ...,
        masses: Optional[np.ndarray] = ...,
    ): ...
    def get_positions(self) -> np.ndarray: ...
    def get_momenta(self) -> np.ndarray: ...
    def get_cell(self) -> np.ndarray: ...
    def get_charges(self) -> np.ndarray: ...
    def get_calculator(self) -> Calculator: ...
    def get_masses(self) -> np.ndarray: ...
    def get_velocities(self) -> np.ndarray: ...
    def get_forces(self) -> np.ndarray: ...
    def get_kinetic_energy(self) -> float: ...
    def get_potential_energy(self) -> float: ...
    def get_initial_charges(self) -> np.ndarray: ...
    def set_positions(self, value: np.ndarray) -> None: ...
    def set_velocities(self, value: np.ndarray) -> None: ...
    def set_momenta(self, value: np.ndarray) -> None: ...
    def set_cell(self, value: np.ndarray) -> None: ...
    def set_calculator(self, value: Calculator) -> None: ...
    def set_pbc(self, value: bool) -> None: ...
    def __len__(self) -> int: ...
    def __iter__(self) -> "Atoms": ...
    def __next__(self) -> Atom: ...
    def append(self, atom: Atom) -> None: ...
    def copy(self) -> Atoms: ...
    @property
    def calc(self) -> Optional[Calculator]: ...
    @calc.setter
    def calc(self, calc: Optional[Calculator]) -> None: ...
