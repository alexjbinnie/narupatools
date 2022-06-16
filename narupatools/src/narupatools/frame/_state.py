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

from typing import Any, Dict, KeysView, Optional, Type, Union

from infinite_sets import InfiniteSet
from narupa.trajectory import FrameData

from narupatools.frame import FrameConverter, FrameKey, get_frame_key
from narupatools.frame._converter import _T


class _StateDataConvert(FrameConverter):
    @classmethod
    def convert_from_frame(
        cls,
        frame: FrameData,
        destination: Union[Type[_T], _T],
        *,
        fields: InfiniteSet[str],
    ) -> _T:
        if destination == StateData:
            state = StateData()
            for key in frame.keys() & fields:  # type: ignore
                state[key] = frame[key]
            return state  # type: ignore
        if isinstance(destination, StateData):
            for key in frame.keys() & fields:  # type: ignore
                destination[key] = frame[key]
            return destination  # type: ignore
        raise NotImplementedError()

    @classmethod
    def convert_to_frame(
        cls, object_: _T, /, *, fields: InfiniteSet[str], existing: Optional[FrameData]
    ) -> FrameData:
        if isinstance(object_, StateData):
            frame = FrameData()
            for key in object_:  # type: ignore
                frame[key] = object_[key]
            return frame
        raise NotImplementedError()


class StateData:
    """Similar to a Narupa frame, but not stored as protobuf."""

    def __init__(self) -> None:
        self._dict: Dict[Union[str, FrameKey], Any] = {}

    def __getitem__(self, key: Union[str, FrameKey]) -> Any:
        if isinstance(key, FrameKey):
            return self._dict[key]
        try:
            return self._dict[get_frame_key(key)]
        except KeyError:
            return self._dict[key]

    def keys(self) -> KeysView[str]:
        """View of keys in the state data."""
        return self._dict.keys()

    def __setitem__(self, key: Union[str, FrameKey], value: Any) -> None:
        if isinstance(key, FrameKey):
            self._dict[key] = key.convert(value)
            return
        try:
            self._dict[key] = get_frame_key(key).convert(value)
        except KeyError:
            self._dict[key] = value

    def __contains__(self, key: Union[str, FrameKey]) -> bool:
        if isinstance(key, FrameKey):
            return key.key in self._dict
        return key in self._dict
