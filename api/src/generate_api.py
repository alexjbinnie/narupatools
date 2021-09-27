from __future__ import annotations

import warnings

import narupatools
from apiparser.info import ModuleInfo, ClassReference, FunctionReference
from apiparser.parse import parse_module
from apiparser.rst import generate_module_rst

with warnings.catch_warnings(record=True) as l:
    analyzed = parse_module(narupatools)

    for w in l:
        print(w.message)

print(analyzed)


modules = {}


def collect_modules(mod: ModuleInfo):
    global modules
    modules[mod.qualified_name] = mod
    for mod in mod.modules:
        collect_modules(mod)


collect_modules(analyzed)


def resolve_refs(mod: ModuleInfo):
    for clss in mod.classes:
        if isinstance(clss, ClassReference):
            clss.target = modules[clss.src_modulename].get_class(clss.src_clssname)
    for func in mod.functions:
        if isinstance(func, FunctionReference):
            func.target = modules[func.value].get_function(func.name)
    for mod in mod.modules:
        resolve_refs(mod)


resolve_refs(analyzed)

generate_module_rst(analyzed, "../../docs/source/api")
