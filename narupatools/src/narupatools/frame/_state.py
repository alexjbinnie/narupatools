from typing import Any, Dict, KeysView, Optional, Type, TypeVar, Union

from infinite_sets import InfiniteSet
from narupa.trajectory import FrameData

from narupatools.override import override

from ._converter import FrameConverter
from .fields import FrameKey, get_frame_key

_T = TypeVar("_T")


class _StateDataConvert(FrameConverter):
    @classmethod
    @override(FrameConverter.convert_from_frame)
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
    @override(FrameConverter.convert_to_frame)
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
