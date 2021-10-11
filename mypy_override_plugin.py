from typing import Optional, Callable, Any

from mypy.nodes import OverloadedFuncDef, Decorator, NameExpr, CallExpr, Block, FuncDef, Argument, Var
from mypy.plugin import Plugin, ClassDefContext

_IGNORED_DECORATORS = [
    "override",
    "serialize_as"
]

def is_ignored(decorator: Any):
    if isinstance(decorator, NameExpr):
        return decorator.name in _IGNORED_DECORATORS
    if isinstance(decorator, CallExpr):
        return decorator.callee.name in _IGNORED_DECORATORS
    return False

class OverridePlugin(Plugin):
    """
    MyPy plugin for ignoring certain decorators on properties.

    Otherwise, MyPy will complain about not permitting decorated properties.
    """

    def clear_decorator(self, decorator):
        decorator.decorators = [dec for dec in decorator.decorators if not is_ignored(dec)]
        decorator.original_decorators = [dec for dec in decorator.original_decorators if
                                         not is_ignored(dec)]

    def modify_context(self, context: ClassDefContext):
        try:
            for expr in context.cls.defs.body:
                if isinstance(expr, Decorator):
                    self.clear_decorator(expr)
                elif isinstance(expr, OverloadedFuncDef):
                    for item in expr.items:
                        if isinstance(item, Decorator):
                            self.clear_decorator(item)
        except Exception as e:
            print(e)

    # use this hook as it is run BEFORE the class is analyzed.
    def get_customize_class_mro_hook(self, fullname: str) -> Optional[Callable[[ClassDefContext], None]]:
        if fullname.startswith("narupatools."):
            return self.modify_context
        if "narupa.trajectory.frame_data.FrameData" == fullname:
            return self.add_monkeypatch

    def add_monkeypatch(self, context: ClassDefContext) -> None:
        block: Block = context.cls.defs
        try:
            block.body.append(FuncDef(name="__getitem__",
                                      arguments=[
                                          Argument(variable=Var("self"), type_annotation=None, initializer=None, kind=0),
                                          Argument(variable=Var("item"), type_annotation=None, initializer=None, kind=0)
                                      ],
                                      body=Block(body=[]),))

            print(block)
        except Exception as e:
            print(e)

def plugin(version: str):
    return OverridePlugin
