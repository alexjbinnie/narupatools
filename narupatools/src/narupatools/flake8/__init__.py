# This file is part of narupatools (https://github.com/alexjbinnie/narupatools).
# Copyright (c) Alex Jamieson-Binnie. All rights reserved.
#
# narupatools is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# narupatools is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with narupatools.  If not, see <http://www.gnu.org/licenses/>.

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
        with open(self.filename) as f:
            first_line = f.readline()
            if not first_line.startswith(
                "# This file is part of narupatools (https://github.com/alexjbinnie/narupatools)."
            ):
                yield 0, 0, "NTL02 Missing copyright notice", type(self)
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
