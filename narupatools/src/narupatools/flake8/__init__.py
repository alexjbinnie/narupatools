"""Narupatools flake8 plugin."""

import ast
from typing import Any, Iterator


class NarupatoolsFlake8Plugin:
    """Flake8 plugin for narupatools."""

    name = "flake8-narupatools"
    version = "1.0.0"

    tree: ast.AST
    filename: str

    def __init__(self, tree: ast.AST, filename: str):
        self.tree = tree
        self.filename = filename

    def run(self) -> Any:  # noqa: D102
        for node in _iter_node(self.tree):
            if isinstance(node, (ast.FunctionDef,)):
                for decorator in node.decorator_list:
                    if (
                        isinstance(decorator, ast.Call)
                        and isinstance(decorator.func, ast.Name)
                        and decorator.func.id == "override"
                    ):
                        pass
                    if isinstance(decorator, ast.Name) and decorator.id == "override":
                        yield decorator.lineno, decorator.col_offset, "NTL01 @override must have one argument", type(
                            self
                        )
            if isinstance(node, ast.ImportFrom) and node.level > 1:
                yield node.lineno, node.col_offset, "NTL02 Cannot import module using ..", type(
                    self
                )


def _iter_node(node: ast.AST) -> Iterator[ast.AST]:
    yield node
    if hasattr(node, "body"):
        for child in node.body:  # type: ignore
            yield from _iter_node(child)
