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

"""Loop which generates frame's from a source at a constant rate."""

from typing import Collection, Protocol

from infinite_sets import InfiniteSet, everything
from narupa.trajectory import FrameData

from narupatools.core._playable import Playable
from narupatools.core.event import Event, EventListener
from narupatools.frame.fields import (
    BondCount,
    BondPairs,
    BoxVectors,
    ChainCount,
    ChainNames,
    KineticEnergy,
    ParticleCount,
    ParticleElements,
    ParticleNames,
    ParticlePositions,
    ParticleResidues,
    ParticleTypes,
    PotentialEnergy,
    ResidueChains,
    ResidueCount,
    ResidueNames,
)
from narupatools.override import override


class ProduceFrameCallback(Protocol):
    """Protocol for a method that can provide a FrameData on demand."""

    def __call__(self, *, fields: InfiniteSet[str]) -> FrameData:
        """
        Produce a new FrameData with the given fields if available.

        :param fields: A set of fields to add to the FrameData if available.
        :return: A FrameData whose only fields are those which were requested and
                 available.
        """
        ...


class OnFrameProducedCallback(Protocol):
    """Callback for when a new frame is produced."""

    def __call__(self, *, frame: FrameData) -> None:
        """
        Called when a new frame is produced.

        :param frame: The FrameData that has been produced.
        """
        ...


DEFAULT_FIELDS = (
    ParticlePositions.key,
    ParticleElements.key,
    ParticleNames.key,
    ParticleTypes.key,
    ParticleResidues.key,
    ResidueNames.key,
    ResidueChains.key,
    ChainNames.key,
    ParticleCount.key,
    ResidueCount.key,
    ChainCount.key,
    BondCount.key,
    BondPairs.key,
    KineticEnergy.key,
    PotentialEnergy.key,
    BoxVectors.key,
)

ALL_FIELDS = everything()


class FrameProducer(Playable):
    """
    Runnable loop that polls a frame source and produces frame at a regular interval.

    A FrameProducer operates on a background thread, and at specified intervals produces
    a FrameData and triggers the on_frame_produced event. However, it only does this if
    the frame has been marked as dirty by using mark_dirty. If mark_dirty is provided a
    set of fields that have changed, then only these fields will then be added to the
    produced frame data.
    """

    _produce: ProduceFrameCallback
    _is_dirty: bool
    _always_dirty: bool
    _fields: InfiniteSet[str]
    _dirty_fields: InfiniteSet[str]
    _on_frame_produced: Event[OnFrameProducedCallback]

    def __init__(
        self,
        produce: ProduceFrameCallback,
        *,
        fields: Collection[str] = DEFAULT_FIELDS,
        frame_interval: float = 1.0 / 30.0,
    ):
        super().__init__(playback_interval=frame_interval)
        self._produce = produce
        self._is_dirty = True
        """Have any fields been marked as dirty since the last frame was produced."""
        self._always_dirty = False
        """Should it be assumed that all fields have changed."""
        self._fields = set(fields)
        """Fields which this producer should read."""
        self._dirty_fields = set(fields)
        """Fields which have been marked as dirty since the last frame was produced."""
        self._on_frame_produced = Event()

    @property
    def always_dirty(self) -> bool:
        """Should it be assumed that all fields are always dirty?"""
        return self._always_dirty

    @always_dirty.setter
    def always_dirty(self, value: bool) -> None:
        self._always_dirty = value

    def add_fields(self, fields: InfiniteSet[str], /) -> None:
        """
        Add one or more fields that should be produced if possible.

        :param fields: Set of fields that should be produced.
        """
        self._fields = self._fields | fields
        self._dirty_fields = self._dirty_fields | fields

    def remove_fields(self, fields: InfiniteSet[str], /) -> None:
        """
        Remove one or more fields that should not be produced.

        :param fields: Set of fields that should no longer be produced.
        """
        self._fields = self._fields - fields
        self._dirty_fields = self._dirty_fields - fields

    @property
    def on_frame_produced(self) -> EventListener[OnFrameProducedCallback]:
        """Event triggered when a new frame is produced."""
        return self._on_frame_produced

    def mark_dirty(self, fields: InfiniteSet[str] = ALL_FIELDS) -> None:
        """
        Indicate that the given fields have been modified.

        :param fields: Set of fields to mark as dirty.
        """
        self._is_dirty = True
        self._dirty_fields = self._dirty_fields | (self._fields & fields)

    @override
    def _advance(self) -> bool:
        if self._is_dirty or self._always_dirty:
            frame = self._produce(fields=self._dirty_fields)
            self._on_frame_produced.invoke(frame=frame)
            self._is_dirty = False
            self._dirty_fields = set()
        return True

    @override
    def _restart(self) -> None:
        self._dirty_fields = everything()
