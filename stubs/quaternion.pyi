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
import numpy as np
import numpy.typing as npt

class quaternion(np.generic):
    def __init__(self, w: float, i: float, j: float, k: float, /): ...
    def conjugate(self) -> quaternion: ...
    def __mul__(self, other: quaternion) -> quaternion: ...
    @property
    def components(self) -> np.ndarray: ...
    x: float
    y: float
    z: float
    w: float
    vec: np.ndarray

one: quaternion

def from_rotation_vector(array: np.ndarray) -> quaternion: ...
def as_rotation_vector(q: quaternion) -> np.ndarray: ...
def as_rotation_matrix(q: quaternion) -> np.ndarray: ...
def as_quat_array(a: np.ndarray) -> npt.NDArray[quaternion]: ...
def from_vector_part(v: npt.NDArray[np.floating]) -> npt.NDArray[quaternion]: ...
