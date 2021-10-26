from narupatools.frame import FrameKey, get_frame_key


class StateData:
    """Similar to a Narupa frame, but not stored as protobuf."""

    def __init__(self):
        self._dict = dict()

    def __getitem__(self, key):
        if isinstance(key, FrameKey):
            return self._dict[key]
        try:
            return self._dict[get_frame_key(key)]
        except KeyError:
            return self._dict[key]

    def __setitem__(self, key, value) -> None:
        if isinstance(key, FrameKey):
            self._dict[key] = key.convert(value)
            return
        try:
            self._dict[key] = get_frame_key(key).convert(value)
        except KeyError:
            self._dict[key] = value

    def __contains__(self, key) -> bool:
        if isinstance(key, FrameKey):
            return key.key in self._dict
        return key in self._dict