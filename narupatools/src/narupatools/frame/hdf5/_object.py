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
