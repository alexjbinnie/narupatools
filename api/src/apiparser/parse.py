import inspect
from types import ModuleType
from typing import Dict, Any, Protocol, Generic

from sphinx.errors import PycodeError
from sphinx.pycode import ModuleAnalyzer
from sphinx.util.inspect import safe_getattr

from .analysis import ModuleAnalysis
from .info import ModuleInfo, ClassReference, ClassInfo, MemberInfo, FunctionReference, FunctionInfo, PropertyInfo
from .util import is_protocol, explore_inheritance, is_class_method, isabstract


def find_instance_members(clss):
    members = {}
    for base_clss in inspect.getmro(clss):
        try:
            modname = base_clss.__module__
            analyzer = ModuleAnalyzer.for_module(modname)
            analyzer.analyze()
        except PycodeError as e:
            continue
        if hasattr(base_clss, "__annotations__"):
            for name, type in base_clss.__annotations__.items():
                if name not in members:
                    members[name] = MemberInfo(clss, name)
                    members[name].is_instancevar = True
                if base_clss != clss and not members[name].inherits:
                    members[name].inherits = [(base_clss, None)]
    return members.values()


def get_module_members(obj: ModuleType) -> Dict[str, Any]:
    members = {}
    for name in dir(obj):
        try:
            members[name] = safe_getattr(obj, name)
        except AttributeError:
            continue
    return members


def parse_method(class_: ClassInfo, name_: str, func_: Any) -> FunctionInfo:
    info = FunctionInfo(name_, func_)
    explore_inheritance(class_.clss, name_, func_, info)
    info.is_class_method = is_class_method(class_, name_, func_)
    info.is_abstract = isabstract(func_)
    return info


def parse_property(class_: ClassInfo, name_: str, property_):
    info = PropertyInfo(name_, property_)
    explore_inheritance(class_.clss, name_, property_, info)
    if info.property.fget:
        sig = inspect.signature(info.property.fget)
        if isinstance(sig.return_annotation, str) and sig.return_annotation.startswith("EventListener"):
            info.is_event = True
    return info


def parse_function(func):
    return FunctionInfo(func.__name__, func)


def parse_class(clss):
    info = ClassInfo(clss)

    # Run the module analyzer provided in Sphinx
    analyzer = ModuleAnalyzer.for_module(clss.__module__)
    analyzer.analyze()

    analysis = ModuleAnalysis(clss.__module__)

    if is_protocol(clss):
        info.is_protocol = True
    for base_class in inspect.getmro(clss):
        if base_class in (clss, object, Protocol):
            continue
        if is_protocol(base_class):
            info.protocols.append(base_class)
        elif base_class is Exception:
            info.is_exception = True
        elif base_class is Generic:
            info.is_generic = True
        else:
            info.base_classes.append(base_class)

    for name, member in inspect.getmembers(clss):
        if inspect.ismethod(member) or inspect.isfunction(member):
            method_info = parse_method(info, name, member)
            if f"{info.name}.{name}" in analyzer.overloads:
                method_info.overloads = analyzer.overloads[f"{info.name}.{name}"]
            info.methods.append(method_info)
        elif name.startswith("__"):
            continue
        elif isinstance(member, property):
            info.properties.append(parse_property(info, name, member))
        else:
            member_info = MemberInfo(info, name)
            member_info.is_classvar = True
            explore_inheritance(clss, name, member, member_info)
            info.members.append(member_info)

    for classname, varname in analyzer.attr_docs:
        if classname != info.name:
            continue
        if any(member.name == varname for member in info.members):
            continue
        if any(member.name == varname for member in info.properties):
            continue
        if any(member.name == varname for member in info.methods):
            continue
        member_info = MemberInfo(info, varname)
        member_info.is_instancevar = True
        member_info.doc = '\n'.join(analyzer.attr_docs[(classname, varname)])
        info.members.append(member_info)

    for classname, varname in analyzer.annotations:
        if classname != info.name:
            continue
        if any(member.name == varname for member in info.members):
            continue
        if any(member.name == varname for member in info.properties):
            continue
        if any(member.name == varname for member in info.methods):
            continue
        member_info = MemberInfo(info, varname)
        member_info.is_instancevar = True
        member_info.annotation = analyzer.annotations[(classname, varname)]
        info.members.append(member_info)

    if info.name in analysis.class_instancevar_assigns:
        for varname in analysis.class_instancevar_assigns[info.name]:
            if any(member.name == varname for member in info.members):
                continue
            if any(member.name == varname for member in info.properties):
                continue
            if any(member.name == varname for member in info.methods):
                continue
            member_info = MemberInfo(info, varname)
            member_info.is_instancevar = True
            info.members.append(member_info)

    return info


def parse_module(module):
    info = ModuleInfo(module)
    analysis = ModuleAnalysis(module.__name__)

    for name, member in get_module_members(module).items():
        # Skip magic methods for now
        # Todo: Implement magic methods so we can tell which are overriden!
        if name.startswith("__"):
            continue
        if inspect.ismodule(member):
            # Skip to avoid recursion
            if member.__name__ == module.__name__:
                continue
            # If submodule
            if member.__name__.startswith(module.__name__):
                info.modules.append(parse_module(member))
            else:
                # Ignore, module imported from elsewhere
                pass
        elif inspect.isclass(member):
            if member.__module__ == module.__name__:
                # Class is defined in this module
                info.classes.append(parse_class(member))
            elif member.__module__.startswith(module.__name__):
                # Class is defined in a submodule
                info.classes.append(ClassReference(member.__module__, name))
            else:
                # For now, ignore classes not defined in a child
                # Todo: Allow class from anywhere
                pass
                # print(f"Skipping class {name} of module {member.__module__}")
        elif inspect.isfunction(member):
            if member.__module__ == module.__name__:
                info.functions.append(parse_function(member))
            elif member.__module__.startswith(module.__name__):
                if name in analysis.imported_symbols:
                    info.functions.append(FunctionReference(info, name, analysis.imported_symbols[name]))
                else:
                    raise ValueError(f"{info.qualified_name} can't work out where {name} came from!")
            else:
                pass
        else:
            if name in analysis.global_assigns:
                member_info = MemberInfo(info, name)
                info.members.append(member_info)
                if ('', name) in analysis.analyzer.attr_docs:
                    member_info.doc = '\n'.join(analysis.analyzer.attr_docs[('', name)])
                if ('', name) in analysis.analyzer.annotations:
                    member_info.annotation = analysis.analyzer.annotations[('', name)]
            elif name in analysis.imported_modules:
                pass
            elif name in analysis.imported_symbols:
                # member_info = MemberReference(info, name, analysis.imported_symbols[name])
                # info.members.append(member_info)
                pass
            else:
                # Member isn't imported, so must have been defined using some black magic
                # Like assigning using setattr. However, it is probably a valid member
                member_info = MemberInfo(info, name)
                info.members.append(member_info)
                if ('', name) in analysis.analyzer.attr_docs:
                    member_info.doc = '\n'.join(analysis.analyzer.attr_docs[('', name)])
    return info
