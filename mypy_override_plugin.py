from typing import Optional, Callable

from mypy.nodes import OverloadedFuncDef, Decorator, NameExpr
from mypy.plugin import Plugin, ClassDefContext


class OverridePlugin(Plugin):
    """
    MyPy plugin for ignoring @override annotations on properties.

    Otherwise, MyPy will complain about not permitting decorated properties.
    """

    def clear_decorator(self, decorator):
        decorator.decorators = [dec for dec in decorator.decorators if
                                not isinstance(dec, NameExpr) or dec.name != "override"]
        decorator.original_decorators = [dec for dec in decorator.original_decorators if
                                         not isinstance(dec, NameExpr) or dec.name != "override"]

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


def plugin(version: str):
    return OverridePlugin
