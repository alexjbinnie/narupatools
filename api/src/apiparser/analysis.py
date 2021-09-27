import ast
import importlib
from typing import List, Dict

from sphinx.pycode import ModuleAnalyzer


class ModuleAnalysis:

    def __init__(self, module):
        self.modname = module
        self.imported_modules: List[str] = []
        """List of modules imported into the global scope."""
        self.imported_symbols: Dict[str, str] = {}
        """Dictionary of objects imported into the global scope, mapping to their source module."""
        self.global_assigns: List[str] = []
        """Symbols assigned to in the global scope."""
        self.class_instancevar_assigns: Dict[str, List[str]] = {}
        """Attributes that are assigned to each class."""

        self._analyzer = ModuleAnalyzer.for_module(module)
        self._analyzer.analyze()

        self.tree = ast.parse(self._analyzer.code)

        self._analyze_assigns(self.tree)
        self._analyze_imports(self.tree)

        self._analyze_instance_var(self.tree)

        self.analyzer = self._analyzer

    def _explore_children(self, node):
        yield node
        if hasattr(node, "body"):
            for child in node.body:
                yield from self._explore_children(child)

    def _analyze_instance_var(self, node):
        if hasattr(node, "body"):
            if isinstance(node, ast.ClassDef):
                classname = node.name
                self.class_instancevar_assigns[classname] = []
                for child in node.body:
                    if isinstance(child, ast.FunctionDef):
                        if not child.args.args or child.args.args[0].arg != "self":
                            continue
                        funcname = child.name
                        for func_node in self._explore_children(child):
                            if isinstance(func_node, ast.Assign):
                                for target in func_node.targets:
                                    if isinstance(target, ast.Attribute):
                                        if isinstance(target.value, ast.Name) and target.value.id == "self":
                                            self.class_instancevar_assigns[classname].append(target.attr)
            else:
                for child in node.body:
                    self._analyze_instance_var(child)

    def _analyze_assigns(self, node):
        if not isinstance(node, (ast.FunctionDef, ast.ClassDef, ast.Lambda)) and hasattr(node, "body"):
            for child in node.body:
                if isinstance(child, (ast.Assign)):
                    for target in child.targets:
                        if isinstance(target, ast.Name):
                            self.global_assigns.append(target.id)
                elif isinstance(child, (ast.AnnAssign)):
                    if isinstance(child.target, ast.Name):
                        self.global_assigns.append(child.target.id)
                self._analyze_assigns(child)

    def _analyze_star_import(self, module_name):
        try:
            mod = importlib.import_module(module_name)
        except ModuleNotFoundError:
            return
        if hasattr(mod, "__all__"):
            for name in mod.__all__:
                self.imported_symbols[name] = module_name
        else:
            for name in dir(mod):
                self.imported_symbols[name] = module_name

    def _analyze_imports(self, node):
        for child in node.body:
            if isinstance(child, ast.ImportFrom):
                mod_name = child.module
                if child.level == 1:
                    mod_name = f"{self.modname}.{mod_name}"
                elif child.level > 1:
                    mod_name = f"{'.'.join(self.modname.split('.')[:-(child.level-1)])}.{mod_name}"
                for name in child.names:
                    if name.name == "*":
                        self._analyze_star_import(mod_name)
                    else:
                        self.imported_symbols[name.name] = mod_name
            elif isinstance(child, ast.Import):
                for name in child.names:
                    self.imported_modules.append(name.name)