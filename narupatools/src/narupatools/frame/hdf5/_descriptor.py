from typing import Any, Optional, Tuple, Type

import numpy as np
import numpy.typing as npt

from narupatools.frame.hdf5._object import _HDF5EditableObject


class HDF5Attribute:
    """Descriptor that stores a value in a HDF5 attribute."""

    def __init__(self, name: str):
        self._name = name

    def __get__(self, instance: _HDF5EditableObject, objtype: Type) -> Any:
        try:
            return instance.hdf5_group._v_attrs[self._name]
        except KeyError as e:
            raise AttributeError(f"Attribute {self._name} is not present.") from e

    def __set__(self, instance: Any, value: Any) -> None:
        if not instance.writable:
            raise ValueError("Trajectory is not writable.")
        instance.hdf5_group._v_attrs[self._name] = value


class HDF5AppendableArray:
    """
    Descriptor representing a HDF5 array that can be appended to.

    When first accessed, the following steps are taken:

    * If the array exists (due to being read from an existing file), it is returned.
    * If
    """

    def __init__(
        self,
        *,
        h5_name: str,
        title: str,
        shape: Tuple[int, ...],
        per_atom: bool = False,
        units: Optional[str],
    ):
        self._h5_name = h5_name
        """Name of the array stored as an HDF5 EArray."""
        self._array_name = f"_{h5_name}_array"
        """Name to store the array as a python attribute."""
        self._title = title
        """Descriptive title of the array."""
        self._per_atom = per_atom
        """Should the size of the array be multiplied by the number of atoms?"""
        self._shape = shape
        """Shape of the numpy array."""
        self._units = units
        """Units of the array."""

    def __get__(self, instance: _HDF5EditableObject, objtype: Type) -> Any:
        try:
            return getattr(instance, self._array_name)
        except AttributeError as e:
            if self._h5_name in instance.hdf5_group:
                value = getattr(instance.hdf5_group, self._h5_name)
            elif instance.writable:
                shape = self._shape
                if self._per_atom:
                    shape = (instance.n_atoms, *shape)
                value = instance.create_earray(
                    name=self._h5_name,
                    title=self._title,
                    shape=shape,
                    units=self._units,
                )
            else:
                raise AttributeError from e
            setattr(instance, self._array_name, value)
            return value

    def as_numpy(self, dtype: npt.DTypeLike = float) -> np.ndarray:
        """Create a property that accesses the array as a numpy array."""

        def get(obj: Any) -> np.ndarray:
            return np.array(self.__get__(obj, type(obj)), dtype=dtype)

        return property(fget=get)  # type: ignore
