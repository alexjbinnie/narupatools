from typing import Any


class AttributeSet:
    def __getitem__(self, name: str) -> Any:
        ...

    def __setitem__(self, key: str, value: Any) -> None:
        ...
