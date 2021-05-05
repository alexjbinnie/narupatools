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

"""Conversion methods for Narupa FrameData."""
from __future__ import annotations

import contextlib
from abc import ABCMeta, abstractmethod
from io import StringIO
from typing import Any, List, Optional, Type, TypeVar, Union, overload
from typing import Dict

import numpy as np
from MDAnalysis.topology.tables import Z2SYMB
from infinite_sets import InfiniteSet, everything
from narupa.trajectory import FrameData

from narupatools.state.typing import Serializable
from .fields import (
    BondPairs,
    ParticleCharges,
    ParticleElements,
    ParticleNames,
    ParticlePositions,
    ParticleResidues,
    ResidueNames,
)
from .frame import NarupaFrame

_T = TypeVar("_T")

_FrameConverters: List[Type[FrameConverter]] = []


class FrameConverter(metaclass=ABCMeta):
    """
    Provides conversions to and from a Narupa FrameData.

    Subclasses are automatically registered to work with the
    :func:`~narupatools.frame.convert` function.
    """

    def __init_subclass__(cls) -> None:
        super().__init_subclass__()
        _FrameConverters.append(cls)

    @classmethod
    @abstractmethod
    def convert_to_frame(  # noqa: D102
        cls, object: _T, *, fields: InfiniteSet[str], existing: Optional[FrameData]
    ) -> FrameData:
        pass

    @classmethod
    @abstractmethod
    def convert_from_frame(  # noqa: D102
        cls, frame: FrameData, type: Union[Type[_T], _T], *, fields: InfiniteSet[str]
    ) -> _T:
        pass


def _wrap_framedata(frame: FrameData) -> NarupaFrame:
    return NarupaFrame(frame.raw)


_TSource = TypeVar("_TSource")
_TTarget = TypeVar("_TTarget")


@overload
def convert(
    source: _TSource, target: Type[_TTarget], *, fields: InfiniteSet[str] = ...
) -> _TTarget:
    ...


@overload
def convert(
    source: _TSource, target: _TTarget, *, fields: InfiniteSet[str] = ...
) -> _TTarget:
    ...


def convert(
    source: _TSource,
    target: Union[Type[_TTarget], _TTarget],
    *,
    fields: InfiniteSet[str] = everything(),
) -> _TTarget:
    """
    Convert between two representations of a molecular system.

    This converts to and from Narupa's FrameData, using it as an intermediary to allow
    conversions between various packages such as ASE and MDAnalysis.

    :param source: Object which contains molecular information.
    :param target: Target type to convert to.
    :param fields: Fields to convert.
    :raises NoConversionDefinedError: No conversion defined between these two objects.
    :return: Object of the requested type.
    """
    target_existing: Optional[Any]
    target_type: Type

    if isinstance(target, type):
        target_existing = None
        target_type = target
    else:
        target_existing = target
        target_type = type(target)

    # Convert from FrameData
    if isinstance(source, FrameData):
        for converter in _FrameConverters:
            with contextlib.suppress(NotImplementedError):
                return converter.convert_from_frame(source, target, fields=fields)
        raise NoConversionDefinedError(source, target)
    # Convert to FrameData
    if target_type in [FrameData, NarupaFrame]:
        for converter in _FrameConverters:
            with contextlib.suppress(NotImplementedError):
                frame_data = converter.convert_to_frame(
                    source, fields=fields, existing=target_existing  # type: ignore
                )
                if target_existing is None and target_type == NarupaFrame:
                    return NarupaFrame(frame_data.raw)  # type: ignore
                else:
                    return frame_data  # type: ignore
        raise NoConversionDefinedError(source, target)
    # Convert in two steps, via Narupa FrameData
    frame: NarupaFrame = convert(source, NarupaFrame)
    return convert(frame, target)


class NoConversionDefinedError(ValueError):
    """Error raised when there is no conversion defined between two objects."""

    def __init__(self, source: Any, target: Any):
        super().__init__(f"No converter exists from {source} to {target}")


class DictConverter(FrameConverter):
    """Converter between a Narupa FrameData and a python dictionary."""

    @classmethod
    def convert_to_frame(
        cls, object: _T, *, fields: InfiniteSet[str], existing: Optional[FrameData]
    ) -> FrameData:  # noqa: D102
        raise NotImplementedError()

    @classmethod
    def convert_from_frame(
        cls, frame: FrameData, type: Union[Type[_T], _T], *, fields: InfiniteSet[str]
    ) -> _T:  # noqa: D102
        if type == Dict:
            return frame_data_to_dictionary(frame, fields=fields)  # type: ignore
        elif isinstance(type, dict):
            return frame_data_to_dictionary(
                frame, fields=fields, existing=type  # type: ignore
            )
        raise NotImplementedError()


def frame_data_to_dictionary(
    frame: FrameData,
    /,
    *,
    fields: InfiniteSet[str] = everything(),
    existing: Optional[Dict[str, Serializable]] = None,
) -> Dict[str, Serializable]:
    """
    Convert a Narupa FrameData to a python dictionary.

    :param frame: Frame to convert to a dictionary.
    :param fields: Fields to copy from the frame.
    :param existing: Existing dictionary to copy into.
    """
    frame_dict: Dict[str, Serializable] = {}
    if existing is not None:
        frame_dict = existing
    for key in frame.values.keys():
        if key in fields:
            frame_dict[key] = frame.values[key]
    for key in frame.arrays.keys():
        if key in fields:
            frame_dict[key] = list(frame.arrays[key])
    return frame_dict


PDB_STD_RESIDUES = frozenset(
    {
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "PYL",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "SEC",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
        "ASX",
        "GLX",
        "UNK",
        "A",
        "C",
        "G",
        "I",
        "U",
        "DA",
        "DC",
        "DG",
        "DI",
        "DT",
        "DU",
        "N",
    }
)
"""
Set of standard residues as defined by the PDB file format.

Any atom in a residue with one of these names should be stored under an ATOM record
instead of a HETATM record.
"""


def frame_to_pdb_string(frame_data: FrameData) -> str:
    """
    Convert a Narupa FrameData to the contents of a PDB file.

    This is useful for interacting with processes which accept a PDB file as input.

    :param frame_data: FrameData to convert.
    :return: Contents of a PDB file as a string.
    """
    file = StringIO()

    positions = ParticlePositions.get(frame_data) * 10.0
    count = len(positions)
    atomnames = ParticleNames.get(frame_data)
    elements = ParticleElements.get(frame_data)
    atomresidues = ParticleResidues.get(frame_data)
    resnames = ResidueNames.get(frame_data)
    try:
        charges = ParticleCharges.get(frame_data)
    except KeyError:
        charges = np.zeros(count)
    file.write("REMARK GENERATED BY NARUPATOOLS\n")

    for index, name, residue, position, element, charge in zip(
        range(count), atomnames, atomresidues, positions, elements, charges
    ):
        record_type = "ATOM  " if resnames[residue] in PDB_STD_RESIDUES else "HETATM"
        file.write(
            f"{record_type}{index + 1:5} {name:4.4} {resnames[residue]:3.3} "
            f"A{residue + 1:4}    {position[0]:8.3f}{position[1]:8.3f}"
            f"{position[2]:8.3f}  1.00  0.00          {Z2SYMB[element]:2.2}"
            f"{int(charge):2}\n"
        )
    bonded: Dict[int, List[int]] = {}
    for bond in BondPairs.get(frame_data):
        if bond[0] not in bonded:
            bonded[bond[0]] = []
        if bond[1] not in bonded:
            bonded[bond[1]] = []
        bonded[bond[0]].append(bond[1])
        bonded[bond[1]].append(bond[0])

    for index in range(count):
        if index not in bonded:
            continue
        line = f"CONECT{index + 1:5}"
        for n in bonded[index]:
            line += f"{n + 1:5}"
        file.write(line + "\n")

    return file.getvalue()
