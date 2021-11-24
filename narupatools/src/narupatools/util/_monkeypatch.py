from typing import Callable, Type


def monkeypatch(baseclass: Type) -> Callable[[Type], None]:
    """Decorate that modifies a base class via monkeypatching."""

    def patch(derivedclass: Type) -> None:
        for key, value in derivedclass.__dict__.items():
            setattr(baseclass, key, value)

    return patch
