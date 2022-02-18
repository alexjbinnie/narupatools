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

"""Describes the base class for anything that supports get_frame."""

from __future__ import annotations

from abc import ABCMeta, abstractmethod
from typing import Any, List, Optional, Protocol, Type

from infinite_sets import InfiniteSet, everything
from narupa.trajectory import FrameData

from narupatools.core.event import EventListener


class OnFieldsChangedCallback(Protocol):
    """Callback for when certain fields of a Frame have changed."""

    def __call__(self, *, fields: InfiniteSet[str]) -> None:  # noqa: D102
        pass


class FrameSource(metaclass=ABCMeta):
    """Base class for object which can create a FrameData."""

    @abstractmethod
    def get_frame(self, *, fields: InfiniteSet[str]) -> FrameData:
        """
        Create a FrameData representation, with the given fields present if available.

        :param fields: Collection of fields that should be added to FrameData if
                       available.
        """


class FrameSourceWithNotify(FrameSource, metaclass=ABCMeta):
    """Base class for objects which inform others when fields are made dirty."""

    @property
    @abstractmethod
    def on_fields_changed(self) -> EventListener[OnFieldsChangedCallback]:  # noqa: D102
        ...


_TRAJ_SOURCE_SUBCLASSES: List[Type[TrajectorySource]] = []


class TrajectorySource(metaclass=ABCMeta):
    """Base class for object which can create multiple FrameData."""

    @abstractmethod
    def get_frame(
        self, *, index: int, fields: InfiniteSet[str] = everything()
    ) -> FrameData:
        """
        Create a FrameData representation, with the given fields present if available.

        :param index: Index of frame to get.
        :param fields: Collection of fields that should be added to FrameData if
                       available.
        """

    @abstractmethod
    def __len__(self) -> int:
        pass

    def __init_subclass__(cls, **kwargs: Any) -> None:
        super().__init_subclass__(**kwargs)
        if cls != TrajectorySource:
            _TRAJ_SOURCE_SUBCLASSES.append(cls)

    @staticmethod
    def create_from_object(obj: Any) -> Optional[TrajectorySource]:
        """Attempt to convert an arbitrary object to a trajectory source."""
        try:
            return TrajectorySource._create_from_object(obj)
        except NotImplementedError:
            return None

    @classmethod
    def _create_from_object(cls, obj: Any) -> TrajectorySource:
        """
        Attempt to convert an arbitrary object to a trajectory source.

        Subclasses can override this to allow automatic conversion of objects to trajectory sources.
        """
        if cls == TrajectorySource:
            for subclass in _TRAJ_SOURCE_SUBCLASSES:
                try:
                    return subclass._create_from_object(obj)
                except NotImplementedError:
                    continue
        raise NotImplementedError
