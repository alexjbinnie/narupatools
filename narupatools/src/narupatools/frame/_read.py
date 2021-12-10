import contextlib
from typing import Any, Union

from narupa.trajectory import FrameData

from narupatools.frame import TrajectorySource

_READER_SUBCLASSES = []


class Reader:

    def load_frame(self, *filenames: str):
        raise ValueError(f"Cannot load {filenames} with {self}.")

    def load_trajectory(self, *filenames: str):
        raise ValueError(f"Cannot load {filenames} with {self}.")

    def __init_subclass__(cls, **kwargs: Any) -> None:
        super().__init_subclass__(**kwargs)  # type: ignore
        if cls != Reader:
            _READER_SUBCLASSES.append(cls)


def load_frame(*filenames: str) -> FrameData:
    for reader in _READER_SUBCLASSES:
        with contextlib.suppress(ValueError):
            return reader.load_frame(*filenames)
    raise ValueError(f"Failed to read files {filenames}")


def load_trajectory(*filenames: str) -> TrajectorySource:
    for reader in _READER_SUBCLASSES:
        with contextlib.suppress(ValueError):
            return reader.load_frame(*filenames)
    raise ValueError(f"Failed to read files {filenames}")


def load(*filenames: str) -> Union[FrameData, TrajectorySource]:
    with contextlib.suppress(ValueError):
        return load_trajectory(*filenames)
    return load_frame(*filenames)