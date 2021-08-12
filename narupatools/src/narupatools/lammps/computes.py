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

"""Defining and extracting computes from LAMMPS."""

from __future__ import annotations

import abc
from typing import Generic, Literal, TypeVar, Union, overload

import numpy as np
import numpy.typing as npt

from narupatools.lammps.constants import VariableDimension, VariableType
from narupatools.lammps.wrapper import LAMMPSWrapper

_TReturnType = TypeVar("_TReturnType")


class ComputeGlobalStyle(Generic[_TReturnType]):
    """Definition of a LAMMPS compute style."""

    def __init__(self, style_id: str, dimension: VariableDimension):
        self.__style_id = style_id
        self.__dimension = dimension

    @overload
    @classmethod
    def define(
        cls, id: str, dimension: Literal[VariableDimension.SCALAR]
    ) -> ComputeGlobalStyle[float]:
        ...

    @overload
    @classmethod
    def define(
        cls,
        id: str,
        dimension: Literal[VariableDimension.VECTOR1D, VariableDimension.ARRAY2D],
    ) -> ComputeGlobalStyle[npt.NDArray[np.float64]]:
        ...

    @overload
    @classmethod
    def define(
        cls, id: str, dimension: VariableDimension
    ) -> Union[ComputeGlobalStyle[npt.NDArray[np.float64]], ComputeGlobalStyle[float]]:
        ...

    @classmethod
    def define(
        cls,
        id: str,
        dimension: VariableDimension,
    ) -> ComputeGlobalStyle:
        """Define a compute style used by LAMMPS."""
        return ComputeGlobalStyle(id, dimension)

    def create(self, lammps: LAMMPSWrapper, id: str) -> ComputeReference[_TReturnType]:
        """Create a new compute using this style."""
        lammps.command(f"compute {id} all {self.__style_id}")
        return ComputeGlobalReference.create(
            lammps,
            id=id,
            dimension=self.__dimension,
        )


KineticEnergy = ComputeGlobalStyle.define(id="ke", dimension=VariableDimension.SCALAR)


class ComputeReference(Generic[_TReturnType], metaclass=abc.ABCMeta):
    """Reference to a previously defined LAMMPS compute."""

    def extract(self) -> _TReturnType:
        """Extract the value of the compute."""
        pass


class ComputeGlobalReference(ComputeReference[_TReturnType]):
    """Reference to a previously defined LAMMPS global compute."""

    def __init__(
        self,
        lammps: LAMMPSWrapper,
        *,
        id: str,
        dimension: VariableDimension,
    ):
        self._lammps = lammps
        self._id = id
        self._dimension = dimension

    @overload
    @classmethod
    def create(
        cls,
        lammps: LAMMPSWrapper,
        *,
        id: str,
        dimension: Literal[VariableDimension.SCALAR],
    ) -> ComputeGlobalReference[float]:
        ...

    @overload
    @classmethod
    def create(
        cls,
        lammps: LAMMPSWrapper,
        *,
        id: str,
        dimension: Literal[VariableDimension.VECTOR1D, VariableDimension.ARRAY2D],
    ) -> ComputeGlobalReference[npt.NDArray[np.float64]]:
        ...

    @overload
    @classmethod
    def create(
        cls,
        lammps: LAMMPSWrapper,
        *,
        id: str,
        dimension: VariableDimension,
    ) -> ComputeGlobalReference:
        ...

    @classmethod
    def create(
        cls,
        lammps: LAMMPSWrapper,
        *,
        id: str,
        dimension: VariableDimension,
    ) -> ComputeGlobalReference:
        """
        Create a compute in the given LAMMPS instance, and return a reference.

        :param lammps: LAMMPS instance to use.
        :param id: The ID to use for this compute.
        :param dimension: The dimension of the data the compute returns.
        :return: A reference to the defined compute, that can be used to extract the value at any time.
        """
        return ComputeGlobalReference(lammps, id=id, dimension=dimension)

    def extract(self) -> _TReturnType:  # noqa: D102
        return self._lammps.extract_global_compute(  # type: ignore[return-value]
            self._id,
            dimension=self._dimension,
        )


class ComputeLocalReference(ComputeReference[_TReturnType]):
    """Reference to a local compute."""

    def __init__(self, lammps: LAMMPSWrapper, *, id: str):
        self._lammps = lammps
        self._id = id

    def extract(self) -> _TReturnType:  # noqa: D102
        return self._lammps.extract_local_compute(self._id)  # type: ignore[return-value]


class ComputeAtomReference(ComputeReference[_TReturnType]):
    """Reference to a per-atom compute, which can be used with `extract_compute` or `gather`."""

    def __init__(
        self,
        lammps: LAMMPSWrapper,
        *,
        id: str,
        type: Literal[VariableType.DOUBLE, VariableType.INTEGER],
        count: int,
    ):
        self._lammps = lammps
        self._id = id
        self._count = count
        self._type = type

    def extract(self) -> _TReturnType:  # noqa: D102
        return self._lammps.extract_atom_compute(self._id)  # type: ignore[return-value]

    def gather(self) -> _TReturnType:
        """Gather this compute across all processors and order by Atom ID."""
        return self._lammps.gather(f"c_{self._id}", self._count)  # type: ignore[return-value]
