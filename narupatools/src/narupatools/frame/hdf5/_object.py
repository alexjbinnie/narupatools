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

from abc import abstractmethod
from typing import Optional, Protocol, Tuple

import tables
from tables import EArray, Float32Atom, Group


class _HDF5EditableObject(Protocol):
    """Define an object backed by a HDF5 group."""

    @property
    @abstractmethod
    def hdf5_group(self) -> Group:
        """HDF5 group to store attributes."""
        pass

    @property
    @abstractmethod
    def writable(self) -> bool:
        """Is this HDF5 object currently writable?"""
        pass

    @property
    @abstractmethod
    def n_atoms(self) -> int:
        """Number of atoms in this object."""
        pass

    def create_earray(
        self,
        *,
        name: str,
        title: str,
        shape: Tuple[int, ...],
        units: Optional[str] = None,
    ) -> EArray:
        """Create a new extendable array."""
        if hasattr(self, "expected_frames"):
            expected_frames = self.expected_frames  # type: ignore
        else:
            expected_frames = 1000
        array = self.hdf5_group._v_file.create_earray(
            self.hdf5_group,
            name=name,
            atom=Float32Atom(),
            shape=(0,) + shape,
            title=title,
            filters=tables.Filters(shuffle=True, complib="zlib", complevel=1),
            expectedrows=expected_frames,
        )
        if units is not None:
            array._v_attrs["units"] = units
        return array  # type: ignore
