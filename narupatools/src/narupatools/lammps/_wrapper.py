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

"""Low-level wrapper around PyLAMMPS to expose some methods in a more pythonic way."""

import abc
import ctypes
import datetime
import typing
from abc import abstractmethod
from ctypes import POINTER, c_char
from typing import Any, Generic, Literal, Optional, Tuple, TypeVar, Union, overload

import numpy as np
import numpy.typing as npt
from lammps import LMP_SIZE_COLS, LMP_SIZE_ROWS, LMP_VAR_ATOM, LMP_VAR_EQUAL, PyLammps

from narupatools.lammps import LAMMPSError
from narupatools.lammps._constants import VariableDimension, VariableStyle, VariableType
from narupatools.lammps._exception_wrapper import catch_lammps_warnings_and_exceptions
from narupatools.lammps.exceptions import (
    ComputeNotFoundError,
    FixNotFoundError,
    GlobalNotFoundError,
    InvalidComputeSpecificationError,
    InvalidThermoKeywordError,
    SettingNotFoundError,
    UnknownAtomPropertyError,
    UnknownPropertyNameError,
    UntestedVersionWarning,
    VariableNotFoundError,
)

SUPPORTED_VERSION = datetime.date(2020, 10, 29)


class LAMMPSWrapper:
    """
    Low level stateless wrapper around LAMMPS.

    This wraps a PyLammps interface (which is already a python wrapper around LAMMPS). This additional
    wrapper has a slightly different interface, with more pythonic behaviour such as raising exceptions
    instead of returning None when invalid keywords are used.
    """

    def __init__(self, pylammps: Optional[PyLammps] = None):
        if pylammps is None:
            pylammps = PyLammps()
        self.__pylammps = pylammps
        self.__lammps = pylammps.lmp
        self.__lammps_clib = pylammps.lmp.lib
        self.__lammps_handle = pylammps.lmp.lmp

        self._check_version()

        self.computes = _Category(self, "compute")
        """Set of currently defined compute IDs."""

        self.dumps = _Category(self, "dump")
        """Set of currently defined dump IDs."""

        self.fixes = _Category(self, "fix")
        """Set of currently defined fix IDs."""

        self.groups = _Category(self, "group")
        """Set of currently defined group IDs."""

        self.molecules = _Category(self, "molecule")
        """Set of currently defined molecule IDs."""

        self.regions = _Category(self, "region")
        """Set of currently defined region IDs."""

        self.variables = _Category(self, "variable")
        """Set of currently defined variable IDs."""

        self.atom_styles = _StyleCategory(self, "atom")
        """Set of possible atom styles."""

        self.integrate_styles = _StyleCategory(self, "integrate")
        """Set of possible integrate styles."""

        self.minimize_styles = _StyleCategory(self, "minimize")
        """Set of possible minimize styles."""

        self.pair_styles = _StyleCategory(self, "pair")
        """Set of possible pair styles."""

        self.bond_styles = _StyleCategory(self, "bond")
        """Set of possible bond styles."""

        self.angle_styles = _StyleCategory(self, "angle")
        """Set of possible angle styles."""

        self.dihedral_styles = _StyleCategory(self, "dihedral")
        """Set of possible dihedral styles."""

        self.improper_styles = _StyleCategory(self, "improper")
        """Set of possible improper styles."""

        self.kspace_styles = _StyleCategory(self, "kspace")
        """Set of possible k-space styles."""

        self.compute_styles = _StyleCategory(self, "compute")
        """Set of possible compute styles."""

        self.fix_styles = _StyleCategory(self, "fix")
        """Set of possible fix styles."""

        self.region_styles = _StyleCategory(self, "region")
        """Set of possible region styles."""

        self.dump_styles = _StyleCategory(self, "dump")
        """Set of possible dump styles."""

        self.command_styles = _StyleCategory(self, "command")
        """Set of possible command styles."""

    @property
    def version_date(self) -> datetime.date:
        """Version of LAMMPS as a :class:`datetime.date`."""
        version = str(self.__lammps.version())
        if len(version) != 8:
            raise ValueError("LAMMPS Version was not of expected format.")
        return datetime.date(int(version[0:4]), int(version[4:6]), int(version[6:8]))

    def _check_version(self) -> None:
        """Check that the current version of LAMMPS is supported."""
        version = self.version_date
        if version != SUPPORTED_VERSION:
            raise UntestedVersionWarning(
                "narupatools was not tested with this version of LAMMPS, and is unsupported."
            )

    def extract_atom_variable(self, name: str) -> npt.NDArray[np.float64]:
        """
        Extract an atom-style variable from the LAMMPS instance.

        An atom-style variable is one defined using a command of the form::

            variable var_name atom ...

        An atom-style variable has one value per atom *local to the current process*. When running with
        MPI, the number of values returned will be *nlocal*. When not using MPI, the length of the array
        will be the total number of atoms.

        .. warning::
           LAMMPS does not differentiate between an equal-style or atom-style variable when extracting.
           If this function is used with an equal-style variable, an array will be returned which does
           not point to the correct data.

        :param name: Name used to define the variable.
        :raises VariableNotFoundError: Variable with the given name does not exist.
        :raises VariableNotFoundError: Variable exists, but is not atom-style or equal-style.
        :return: NumPy float array containing the value of the variable for each atom on the current process.
        """
        if not isinstance(name, str):
            raise VariableNotFoundError(name)
        value = self.__lammps.numpy.extract_variable(name, vartype=LMP_VAR_ATOM)
        if value is None or value.shape == ():  # type: ignore[union-attr]
            raise VariableNotFoundError(name)
        return value  # type: ignore[return-value]

    def extract_variable(self, name: str) -> float:
        """
        Extract an equal-style variable from the LAMMPS instance.

        An equal-style variable is one defined using a command of the form::

            variable var_name equal ...

        An equal-style variable has a single float value.

        .. warning::
           LAMMPS does not differentiate between an equal-style or atom-style variable when extracting.
           If this function is used with an atom-style variable, the value may be incorrect.

        :param name: Name used to define the variable.
        :raises VariableNotFoundError: Variable with the given name does not exist.
        :raises VariableNotFoundError: Variable exists, but is not atom-style or equal-style.
        :return: Current value of the variable as a float.
        """
        if not isinstance(name, str):
            raise VariableNotFoundError(name)
        value = self.__lammps.numpy.extract_variable(name, vartype=LMP_VAR_EQUAL)
        if value is None:
            raise VariableNotFoundError(name)
        return value  # type: ignore[return-value]

    def extract_setting(self, name: str) -> int:
        """
        Retrieve or compute a global setting in the LAMMPS instance.

        A list of valid settings are hardcoded in the LAMMPS library.

        :param name: Name of the global setting.
        :raises SettingNotFoundError: Setting with the given name does not exist.
        :return: Value of the global setting.
        """
        value = self.__lammps.extract_setting(name)
        if value == -1:
            raise SettingNotFoundError(name)
        return value

    def extract_global(self, name: str) -> Union[str, float, int]:
        """
        Extract a global property defined in the LAMMPS instance.

        Properties can be either an integer, a float or a string. The type is queried
        directly from LAMMPS using the `extract_global_datatype` function.

        :param name: Name of the global property.
        :raises GlobalNotFoundError: Global with the given name was not found.
        :return: Value of the global property, with the correct python type.
        """
        if not isinstance(name, str):
            raise GlobalNotFoundError(name)

        type = self.__lammps.extract_global_datatype(name)

        if type == -1:
            raise GlobalNotFoundError(name)

        if type == VariableType.STRING:
            self.__lammps_clib.lammps_extract_global.restype = POINTER(c_char)
            ptr = self.__lammps_clib.lammps_extract_global(
                self.__lammps_handle, name.encode()
            )
            return ctypes.string_at(ptr).decode("ascii")
        else:
            value = self.__lammps.extract_global(name, type)
            if value is None:
                raise GlobalNotFoundError(name)
            return value

    @overload
    def extract_global_compute(
        self,
        compute_id: str,
        *,
        dimension: Literal[VariableDimension.SCALAR],
    ) -> float:
        ...

    @overload
    def extract_global_compute(
        self,
        compute_id: str,
        *,
        dimension: Literal[VariableDimension.VECTOR1D, VariableDimension.ARRAY2D],
    ) -> npt.NDArray[np.float64]:
        ...

    @overload
    def extract_global_compute(
        self,
        compute_id: str,
        *,
        dimension: VariableDimension,
    ) -> Union[float, npt.NDArray[np.float64]]:
        ...

    def extract_global_compute(
        self,
        compute_id: str,
        *,
        dimension: VariableDimension,
    ) -> Union[float, npt.NDArray[np.float64]]:
        """
        Extract value of a compute which returns a global variable.

        As computes can return one or more of a scalar, vector (1D) or array (2D) value, the
        dimensionality of the data must also be specified.

        :param compute_id: ID of the compute.
        :param dimension: Dimension of the requested return value.
        :raises ComputeNotFoundError: Key does not match any known computes.
        :raises InvalidComputeSpecificationError: Compute exists but isn't global, or does not
                                                  match dimension provided.
        :return: Value of the compute.
        """
        with catch_lammps_warnings_and_exceptions():
            try:
                value = self.__lammps.numpy.extract_compute(
                    compute_id, VariableStyle.GLOBAL, dimension
                )
            except ValueError as e:
                if e.args[0] == "NULL pointer access":
                    if compute_id not in self.computes:
                        raise ComputeNotFoundError(compute_id)
                    else:
                        raise InvalidComputeSpecificationError(
                            compute_id, VariableStyle.GLOBAL, dimension
                        )
                else:
                    raise e
        if value is None:
            if compute_id not in self.computes:
                raise ComputeNotFoundError(compute_id)
            else:
                raise InvalidComputeSpecificationError(
                    compute_id, VariableStyle.GLOBAL, dimension
                )
        if dimension == VariableDimension.VECTOR1D:
            return value.flatten()  # type: ignore[union-attr]
        return value

    def extract_local_compute(self, compute_id: str) -> npt.NDArray[np.float64]:
        """
        Extract value of a compute which returns a local variable.

        :param compute_id: ID of the compute.
        :raises ComputeNotFoundError: Key does not match any known computes.
        :return: Value of the compute.
        """
        try:
            nrows = self.__lammps.extract_compute(
                compute_id, VariableStyle.LOCAL, LMP_SIZE_ROWS
            )
            ncols = self.__lammps.extract_compute(
                compute_id, VariableStyle.LOCAL, LMP_SIZE_COLS
            )
        except ValueError as e:
            if e.args[0] == "NULL pointer access":
                raise ComputeNotFoundError(compute_id)
            else:
                raise e

        if ncols == 0:
            type: Literal[
                VariableDimension.VECTOR1D, VariableDimension.ARRAY2D
            ] = VariableDimension.VECTOR1D
        else:
            type = VariableDimension.ARRAY2D

        value = self.__lammps.extract_compute(compute_id, VariableStyle.LOCAL, type)

        if type == VariableDimension.VECTOR1D:
            return self._to_numpy(value, (nrows,))
        elif type == VariableDimension.ARRAY2D:
            return self._to_numpy(value, (nrows, ncols))

    def extract_atom_compute(self, compute_id: str) -> npt.NDArray[np.float64]:
        """
        Extract value of a compute which returns a per-atom variable.

        :param compute_id: ID of the compute.
        :raises ComputeNotFoundError: Key does not match any known computes.
        :return: Value of the compute.
        """
        try:
            ncols = self.__lammps.extract_compute(
                compute_id, VariableStyle.ATOM, LMP_SIZE_COLS
            )
        except ValueError as e:
            if e.args[0] == "NULL pointer access":
                raise ComputeNotFoundError(compute_id)
            else:
                raise e

        natoms = self.__lammps.get_natoms()

        if ncols == 0:
            type: Literal[
                VariableDimension.VECTOR1D, VariableDimension.ARRAY2D
            ] = VariableDimension.VECTOR1D
        else:
            type = VariableDimension.ARRAY2D

        value = self.__lammps.extract_compute(compute_id, VariableStyle.ATOM, type)

        if type == VariableDimension.VECTOR1D:
            return self._to_numpy(value, (natoms,))
        elif type == VariableDimension.ARRAY2D:
            return self._to_numpy(value, (natoms, ncols))

    def _to_numpy(self, raw_ptr: Any, shape: Tuple[int, ...]) -> np.ndarray:
        if len(shape) == 1:
            ptr = ctypes.cast(raw_ptr, POINTER(ctypes.c_double * shape[0]))
        else:
            ptr = ctypes.cast(
                raw_ptr[0], POINTER(ctypes.c_double * (shape[0] * shape[1]))
            )

        array = np.frombuffer(ptr.contents)
        array.shape = shape
        array.flags.writeable = False
        return array  # type: ignore

    def extract_atom_fix(self, fix_id: str) -> npt.NDArray[np.float64]:
        """
        Extract value of a fix which returns a per-atom variable.

        :param fix_id: ID of the fix.
        :raises ComputeNotFoundError: Key does not match any known fixes.
        :return: Value of the fix.
        """
        try:
            ncols = self.__lammps.extract_fix(fix_id, VariableStyle.ATOM, LMP_SIZE_COLS)
        except ValueError as e:
            if e.args[0] == "NULL pointer access":
                raise FixNotFoundError(fix_id)
            else:
                raise e

        natoms = self.__lammps.get_natoms()

        if ncols == 0:
            type: Literal[
                VariableDimension.VECTOR1D, VariableDimension.ARRAY2D
            ] = VariableDimension.VECTOR1D
        else:
            type = VariableDimension.ARRAY2D

        value = self.__lammps.extract_fix(fix_id, VariableStyle.ATOM, type)

        if type == VariableDimension.VECTOR1D:
            return self._to_numpy(value, (natoms,))
        elif type == VariableDimension.ARRAY2D:
            return self._to_numpy(value, (natoms, ncols))

    def gather_atoms(
        self,
        key: str,
        dimension: int,
    ) -> Union[npt.NDArray[np.float64], npt.NDArray[np.int32]]:
        """
        Gather a per-atom property from all the processors.

        :param key: Id of the property.
        :param dimension: Dimension of the variable.
        :raises UnknownAtomPropertyError: Property was not found in the simulation.
        :return: Value of the atom property.
        """
        type = self.__lammps.extract_atom_datatype(key)
        if type == -1:
            raise UnknownAtomPropertyError(key)
        try:
            if type in [VariableType.DOUBLE, VariableType.DOUBLE_ARRAY]:
                dtype: npt.DTypeLike = np.float64
                ctype: Any = ctypes.c_double
                type = VariableType.DOUBLE
            elif type in [VariableType.INTEGER, VariableType.INTEGER_ARRAY]:
                dtype = np.int32
                ctype = ctypes.c_int
                type = VariableType.INTEGER
            else:
                raise LAMMPSError(f"Invalid variable type {type}")
            with catch_lammps_warnings_and_exceptions():
                natoms = self.__lammps.get_natoms()
                data = ((dimension * natoms) * ctype)()

                self.__lammps_clib.lammps_gather_atoms(
                    self.__lammps_handle, key.encode(), type, dimension, data
                )

            if dimension == 1:
                array = np.array(data, dtype=dtype).reshape((-1,))
            else:
                array = np.array(data, dtype=dtype).reshape((-1, dimension))
            array.flags.writeable = False
            return array
        except UnknownPropertyNameError:
            # Reraise as a different error so we know what key caused the error.
            raise UnknownAtomPropertyError(key)

    def gather_compute(self, compute_id: str) -> npt.NDArray[np.float64]:
        """
        Gather the value of a per-atom compute across all processors, sorted by atom ID.

        :param compute_id: ID of the compute to gather.
        :return: Value of the compute for each atom, ordered by atom ID.
        """
        try:
            ncols = self.__lammps.extract_compute(
                compute_id, VariableStyle.ATOM, LMP_SIZE_COLS
            )
        except ValueError as e:
            if e.args[0] == "NULL pointer access":
                raise ComputeNotFoundError(compute_id)
            else:
                raise e
        return self.gather("c_" + compute_id, ncols)

    def gather_fix(self, fix_id: str) -> npt.NDArray[np.float64]:
        """
        Gather the value of a per-atom fix across all processors, sorted by atom ID.

        :param fix_id: ID of the fix to gather.
        :return: Value of the fix for each atom, ordered by atom ID.
        """
        try:
            ncols = self.__lammps.extract_fix(fix_id, VariableStyle.ATOM, LMP_SIZE_COLS)
        except ValueError as e:
            if e.args[0] == "NULL pointer access":
                raise FixNotFoundError(fix_id)
            else:
                raise e
        return self.gather("f_" + fix_id, ncols)

    def gather(
        self,
        key: str,
        dimension: int,
    ) -> np.ndarray:
        """See :meth:`lammps.gather`."""
        if dimension == 0:
            dimension = 1

        if key[:2] in ["d_", "c_", "f_"]:
            type = VariableType.DOUBLE
        elif key[:2] in ["i_"]:
            type = VariableType.INTEGER
        else:
            raw = self.__lammps.extract_atom_datatype(key)
            if raw == -1:
                raise UnknownAtomPropertyError(key)
            type = VariableType(raw)

        if type in [VariableType.DOUBLE, VariableType.DOUBLE_ARRAY]:
            dtype: npt.DTypeLike = np.float64
            ctype: Any = ctypes.c_double
        elif type in [VariableType.INTEGER, VariableType.INTEGER_ARRAY]:
            dtype = np.int32
            ctype = ctypes.c_int
        else:
            raise LAMMPSError(f"Invalid variable type {type}")

        with catch_lammps_warnings_and_exceptions():
            natoms = self.__lammps.get_natoms()
            data = ((dimension * natoms) * ctype)()

            self.__lammps_clib.lammps_gather(
                self.__lammps_handle, key.encode(), type, dimension, data
            )

        if dimension == 1:
            array = np.array(data, dtype=dtype).reshape((-1,))
        else:
            array = np.array(data, dtype=dtype).reshape((-1, dimension))
        array.flags.writeable = False
        return array

    def scatter_atoms(
        self, name: str, type: VariableType, dimensions: int, value: np.ndarray
    ) -> None:
        """See :meth:`lammps.scatter_atoms`."""
        n_atoms = self.__lammps.get_natoms()
        with catch_lammps_warnings_and_exceptions():
            self.__lammps.scatter_atoms(
                name, type, dimensions, _to_ctypes(value, n_atoms)
            )

    def get_thermo(self, keyword: str) -> float:
        """
        Evaluate a thermo keyword as defined in LAMMPS.

        This returns 0.0 if the keyword is invalid.

        :param keyword: Thermo keyword to compute.
        :raises InvalidThermoKeywordError: Thermo keyword is not a string.
        :return: Value of the keyword, or 0.0 if the keyword is unrecognized.
        """
        value = self.__lammps.get_thermo(keyword)
        if value is None:
            raise InvalidThermoKeywordError(keyword)
        return value

    def close(self) -> None:
        """Close the LAMMPS instance."""
        self.__pylammps.close()

    def file(self, filename: str) -> None:
        """See :meth:`PyLammps.file`."""
        self.__pylammps.file(filename)

    def _has_id(self, category: str, id: str) -> bool:
        return self.__lammps.has_id(category, id)

    def _has_style(self, category: str, id: str) -> bool:
        return self.__lammps.has_style(category, id)

    def _id_count(self, category: str) -> int:
        return self.__lammps_clib.lammps_id_count(  # type:ignore[no-any-return]
            self.__lammps_handle, category.encode()
        )

    def _style_count(self, category: str) -> int:
        return self.__lammps_clib.lammps_style_count(  # type:ignore[no-any-return]
            self.__lammps_handle, category.encode()
        )

    def _id_name(self, category: str, index: int) -> str:
        sb = ctypes.create_string_buffer(100)
        result = self.__lammps_clib.lammps_id_name(
            self.__lammps_handle, category.encode(), index, sb, 100
        )
        if result == 0:
            raise IndexError(index)
        return sb.value.decode()

    def _style_name(self, category: str, index: int) -> str:
        sb = ctypes.create_string_buffer(100)
        result = self.__lammps_clib.lammps_style_name(
            self.__lammps_handle, category.encode(), index, sb, 100
        )
        if result == 0:
            raise IndexError(index)
        return sb.value.decode()

    def command(self, command: str) -> None:
        """
        Run an arbitrary LAMMPS command.

        Running it through this interface means that all previously cached data such as
        positions etc. are cleared. If you are sure running the command will not affect
        any data, set clear_cache to False.

        :param command: LAMMPS command to run.
        :param clear_cache: Should all cached information be cleared.
        """
        with catch_lammps_warnings_and_exceptions():
            self.__pylammps.command(command)


class _StyleCategory(typing.Sequence[str]):
    def __init__(self, lammps: LAMMPSWrapper, category: str):
        self.__lammps = lammps
        self.__category = category

    def __contains__(self, id: object) -> bool:
        return isinstance(id, str) and self.__lammps._has_style(self.__category, id)

    def __len__(self) -> int:
        return self.__lammps._style_count(self.__category)

    def __getitem__(self, index: object) -> str:
        if not isinstance(index, int):
            raise IndexError(index)
        return self.__lammps._style_name(self.__category, index)


class _Category(typing.Sequence[str]):
    def __init__(self, lammps: LAMMPSWrapper, category: str):
        self.__lammps = lammps
        self.__category = category

    def __contains__(self, id: object) -> bool:
        return isinstance(id, str) and self.__lammps._has_id(self.__category, id)

    def __len__(self) -> int:
        return self.__lammps._id_count(self.__category)

    def __getitem__(self, index: object) -> str:
        if not isinstance(index, int):
            raise IndexError(index)
        return self.__lammps._id_name(self.__category, index)

    def generate_id(self) -> str:
        x = 1
        ids = list(self)
        while str(x) in ids:
            x += 1
        return str(x)


_TReturnType = TypeVar("_TReturnType")


class Extractable(Generic[_TReturnType], metaclass=abc.ABCMeta):
    """Define something which can be extracted from LAMMPS."""

    @abstractmethod
    def extract(self, lammps: LAMMPSWrapper) -> _TReturnType:
        """Extract the value from a LAMMPS instance."""
        pass


def _to_ctypes(array: np.ndarray, natoms: int) -> ctypes.Array:
    n3 = 3 * natoms
    x = (n3 * ctypes.c_double)()
    for i, f in enumerate(array.flat):
        x[i] = f
    return x
