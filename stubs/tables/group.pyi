from typing import Any

class Group:
    def __getattr__(self, name: str) -> Any: ...
    _v_nchildren: int

class RootGroup(Group): ...
