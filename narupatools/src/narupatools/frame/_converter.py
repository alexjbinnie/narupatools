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
from typing import Any, Dict, List, Optional, Type, TypeVar, Union, overload

from infinite_sets import InfiniteSet, everything
from narupa.trajectory import FrameData

from narupatools.state.typing import Serializable

from ..override import override

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
        cls, object_: _T, /, *, fields: InfiniteSet[str], existing: Optional[FrameData]
    ) -> FrameData:
        pass

    @classmethod
    @abstractmethod
    def convert_from_frame(  # noqa: D102
        cls,
        frame: FrameData,
        destination: Union[Type[_T], _T],
        *,
        fields: InfiniteSet[str],
    ) -> _T:
        pass


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
    if target_type is FrameData:
        for converter in _FrameConverters:
            with contextlib.suppress(NotImplementedError):
                return converter.convert_to_frame(
                    source, fields=fields, existing=target_existing  # type: ignore
                )
        raise NoConversionDefinedError(source, target)
    # Convert in two steps, via Narupa FrameData
    frame: FrameData = convert(source, FrameData)
    return convert(frame, target)


class NoConversionDefinedError(ValueError):
    """Error raised when there is no conversion defined between two objects."""

    def __init__(self, source: Any, target: Any):
        super().__init__(f"No converter exists from {source} to {target}")


class DictConverter(FrameConverter):
    """Converter between a Narupa FrameData and a python dictionary."""

    @classmethod
    @override(FrameConverter.convert_to_frame)
    def convert_to_frame(  # noqa: D102
        cls, object_: _T, /, *, fields: InfiniteSet[str], existing: Optional[FrameData]
    ) -> FrameData:
        raise NotImplementedError

    @classmethod
    @override(FrameConverter.convert_from_frame)
    def convert_from_frame(  # noqa: D102
        cls,
        frame: FrameData,
        destination: Union[Type[_T], _T],
        *,
        fields: InfiniteSet[str],
    ) -> _T:
        if destination == Dict:
            return frame_data_to_dictionary(frame, fields=fields)  # type: ignore
        elif isinstance(destination, dict):
            return frame_data_to_dictionary(
                frame, fields=fields, existing=destination  # type: ignore
            )
        raise NotImplementedError


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
