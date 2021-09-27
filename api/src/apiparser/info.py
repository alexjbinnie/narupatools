from __future__ import annotations

import inspect
from types import ModuleType
from typing import Union, List, TypeVar, Generic, Optional, Any

from sphinx.pycode import ModuleAnalyzer

from .util import prepend_tab


class MemberInfo:
    def __init__(self, parent: Union[ModuleInfo, ClassInfo], name: str):
        self.parent = parent
        if inspect.ismodule(parent):
            analyzer = ModuleAnalyzer.for_file(parent.module.__file__, parent.module.__name__)
            analyzer.analyze()
            try:
                self.doc = "".join(analyzer.attr_docs[('', name)])
            except KeyError:
                self.doc = None
        elif inspect.isclass(parent):
            analyzer = ModuleAnalyzer.for_module(parent.clss.__module__)
            analyzer.analyze()
            try:
                self.doc = "".join(analyzer.attr_docs[(parent.clss.__name__, name)])
            except KeyError:
                self.doc = None
        else:
            self.doc = None
        self.annotation = None
        self.name = name
        self.is_classvar = False
        self.is_instancevar = False
        self.implements = []
        self.overrides = []
        self.inherits = []

    @property
    def qualified_name(self):
        return f"{self.parent.qualified_name}.{self.name}"

    def __str__(self):
        s = f"Member {self.name}"
        if self.doc:
            s += " <documented>"
        if self.annotation:
            s += " <annotated>"
        if self.is_classvar:
            s += " <classvar>"
        if self.is_instancevar:
            s += " <instancevar>"
        if self.implements:
            s += f" <implements from {self.implements[0][0].__name__}>"
        if self.overrides:
            s += f" <overrides from {self.overrides[0][0].__name__}>"
        if self.inherits:
            s += f" <inherits from {self.inherits[0][0].__name__}>"
        return s


class ModuleInfo:
    def __init__(self, module: ModuleType):
        self.qualified_name = module.__name__
        self.name = self.qualified_name.split(".")[-1]
        self.module = module
        try:
            self.doc = module.__doc__
        except AttributeError:
            self.doc = None
        self.modules: List[ModuleInfo] = []
        self.classes: List[Union[ClassInfo, ClassReference]] = []
        self.functions: List[Union[FunctionInfo, FunctionReference]] = []
        self.members: List[Union[MemberInfo, MemberReference]] = []

    def get_class(self, name: str) -> Union[ClassInfo, ClassReference]:
        for clss in self.classes:
            if clss.name == name:
                return clss
        raise KeyError

    def get_function(self, name: str) -> Union[FunctionInfo, FunctionReference]:
        for func in self.functions:
            if func.name == name:
                return func
        raise KeyError(f"Failed to find function {name} in module {self.name}")

    def __str__(self):
        s = f"Module {self.name}"
        if self.doc:
            s += " <documented>"
        for module in self.modules:
            s += prepend_tab("\n" + str(module))
        for clss in self.classes:
            s += prepend_tab("\n" + str(clss))
        for func in self.functions:
            s += prepend_tab("\n" + str(func))
        for member in self.members:
            s += prepend_tab("\n" + str(member))
        return s


class ClassInfo:
    def __init__(self, clss):
        self.name = clss.__name__
        self.qualified_name = f"{clss.__module__}.{clss.__name__}"
        self.clss = clss
        self.is_protocol = False
        self.is_generic = False
        self.is_exception = False
        try:
            self.doc = clss.__doc__
        except AttributeError:
            self.doc = None
        self.protocols = []
        self.base_classes = []
        self.methods = []
        self.members: List[MemberInfo] = []
        self.properties: List[PropertyInfo] = []

    def has_method(self, name: str):
        return any(method.name == name for method in self.methods)

    def __str__(self):
        s = f"Class {self.name}"
        if self.doc:
            s += " <documented>"
        if self.is_protocol:
            s += " <protocol>"
        if self.is_generic:
            s += " <generic>"
        if self.is_exception:
            s += " <exception>"
        for protocol in self.protocols:
            s += prepend_tab(f"\nImplements {protocol.__name__}")
        for base_class in self.base_classes:
            s += prepend_tab(f"\nInherits {base_class.__name__}")
        for property_ in self.properties:
            s += prepend_tab("\n" + str(property_))
        for func in self.methods:
            s += prepend_tab("\n" + str(func))
        for func in self.members:
            s += prepend_tab("\n" + str(func))
        return s


class FunctionInfo:

    def __init__(self, name, method):
        self.name = name
        self.method = method
        self.is_class_method = False
        self.is_abstract = False
        try:
            self.doc = method.__doc__
        except AttributeError:
            self.doc = None
        self.implements = []
        self.overrides = []
        self.inherits = []
        self.overloads = []

    @property
    def signature(self):
        return inspect.signature(self.method)

    def __str__(self):
        s = f"Method {self.name}"
        if self.doc:
            s += " <documented>"
        if self.is_class_method:
            s += " <classmethod>"
        if self.is_abstract:
            s += " <abstract>"
        if self.overloads:
            s += " <overloaded>"
        if self.implements:
            s += f" <implements from {self.implements[0][0].__name__}>"
        if self.overrides:
            s += f" <overrides from {self.overrides[0][0].__name__}>"
        if self.inherits:
            s += f" <inherits from {self.inherits[0][0].__name__}>"

        return s


class PropertyInfo:
    def __init__(self, name, property_):
        self.name = name
        self.property = property_
        self.is_abstract = False
        self.is_event = False

        self.implements = []
        self.overrides = []
        self.inherits = []

    @property
    def doc(self):
        if self.property.__doc__ is not None:
            return self.property.__doc__
        for _, implement in self.implements:
            if implement.__doc__ is not None:
                return implement.__doc__
        for _, implement in self.overrides:
            if implement.__doc__ is not None:
                return implement.__doc__
        for _, implement in self.inherits:
            if implement.__doc__ is not None:
                return implement.__doc__
        return None

    def __str__(self):
        s = f"Property {self.name}"
        if self.doc:
            s += " <documented>"
        if self.property.fget:
            s += " <getter>"
        if self.property.fset:
            s += " <setter>"
        if self.is_abstract:
            s += " <abstract>"
        if self.implements:
            s += f" <implements from {self.implements[0][0].__name__}>"
        if self.overrides:
            s += f" <overrides from {self.overrides[0][0].__name__}>"
        if self.inherits:
            s += f" <inherits from {self.inherits[0][0].__name__}>"
        return s


_TInfo = TypeVar("_TInfo", bound=Union[MemberInfo, ClassInfo, FunctionInfo])


class Reference(Generic[_TInfo]):
    def __init__(self):
        self.target: Optional[_TInfo] = None

    def resolve(self, target: _TInfo):
        self.target = target

    def __getattr__(self, item):
        if self.target is None:
            raise AttributeError(f"Reference not resolved ({self}) when getting {item}")
        return getattr(self.target, item)


class ClassReference(Reference[FunctionInfo]):
    def __init__(self, src_modulename: str, src_clssname:str):
        super().__init__()
        self.src_modulename = src_modulename
        self.src_clssname = src_clssname

    def __str__(self):
        if self.target is not None:
            return f"ClassReference to {self.target}"
        else:
            return f"ClassReference to {self.src_clssname} in {self.src_modulename}"


class MemberReference(Reference[FunctionInfo]):

    def __init__(self, referee, name, srcmodule):
        super().__init__()
        self.referee = referee
        self.name = name
        self.srcmodule = srcmodule
        self.target = None

    def __str__(self):
        if self.target is not None:
            return f"MemberReference to {self.target}"
        else:
            return f"MemberReference in {self.referee.name} to {self.name} from {self.srcmodule}"


class FunctionReference(Reference[FunctionInfo]):

    def __init__(self, referee_module: ModuleInfo, func_name: str, value: Any):
        super().__init__()
        self.referee_module = referee_module
        self.name = func_name
        self.value = value

    def __str__(self):
        if self.target is not None:
            return f"FunctionReference to {self.target}"
        else:
            return f"FunctionReference in {self.referee_module.name} to {self.name} from {self.value}"
