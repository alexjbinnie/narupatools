"""Plugin for MyPy."""

from typing import Any, Callable, Optional

from mypy.nodes import (
    Argument,
    Block,
    CallExpr,
    Decorator,
    FuncDef,
    NameExpr,
    OverloadedFuncDef,
    Var,
)
from mypy.plugin import ClassDefContext, Plugin

_IGNORED_DECORATORS = ["override", "serialize_as"]


def is_ignored(decorator: Any):
    """Should the given decorator be ignored by MyPy."""
    if isinstance(decorator, NameExpr):
        return decorator.name in _IGNORED_DECORATORS
    if isinstance(decorator, CallExpr):
        return decorator.callee.name in _IGNORED_DECORATORS
    return False


class NarupatoolsPlugin(Plugin):
    """
    MyPy plugin for narupatools.

    This allows certain decorators on properties to be ignored, to avoid MyPy complaining.

    It also tricks MyPy into allowing __getitem__ and __setitem__ for FrameData, as these methods are
    monkeypatched by narupatools.
    """

    def remove_ignored_decorators(self, decorator):
        """Remove decorators which should be ignored."""
        decorator.decorators = [
            dec for dec in decorator.decorators if not is_ignored(dec)
        ]
        decorator.original_decorators = [
            dec for dec in decorator.original_decorators if not is_ignored(dec)
        ]

    def remove_decorators(self, context: ClassDefContext):
        """Callback that removes ignored decorators from narupatools classes."""
        for expr in context.cls.defs.body:
            if isinstance(expr, Decorator):
                self.remove_ignored_decorators(expr)
            elif isinstance(expr, OverloadedFuncDef):
                for item in expr.items:
                    if isinstance(item, Decorator):
                        self.remove_ignored_decorators(item)

    def get_customize_class_mro_hook(  # noqa: D102
            self, fullname: str
    ) -> Optional[Callable[[ClassDefContext], None]]:
        if fullname.startswith("narupatools."):
            return self.remove_decorators
        if fullname == "narupa.trajectory.frame_data.FrameData":
            return self.add_monkeypatched_framedata

    def add_monkeypatched_framedata(self, context: ClassDefContext) -> None:
        """Modify FrameData definition to include monkey-patched functions."""
        setter = FuncDef(
            name="__setitem__",
            arguments=[
                Argument(
                    variable=Var("self"), type_annotation=None, initializer=None, kind=0
                ),
                Argument(
                    variable=Var("key"), type_annotation=None, initializer=None, kind=0
                ),
                Argument(
                    variable=Var("value"),
                    type_annotation=None,
                    initializer=None,
                    kind=0,
                ),
            ],
            body=Block(body=[]),
        )

        context.cls.defs.body.append(setter)

        getter = FuncDef(
            name="__getitem__",
            arguments=[
                Argument(
                    variable=Var("self"), type_annotation=None, initializer=None, kind=0
                ),
                Argument(
                    variable=Var("key"), type_annotation=None, initializer=None, kind=0
                ),
            ],
            body=Block(body=[]),
        )

        context.cls.defs.body.append(getter)


def plugin(version: str):
    """Entry point of the mypy plugin."""
    return NarupatoolsPlugin
