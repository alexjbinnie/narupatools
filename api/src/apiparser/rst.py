import inspect
import os
from contextlib import redirect_stdout
from typing import Union, get_origin, get_args, Dict

from .info import ClassInfo, ClassReference, MemberReference, MemberInfo, FunctionReference, FunctionInfo, \
    PropertyInfo, ModuleInfo
from .util import prepend_spaces, indent


def generate_module_member_rst(member: Union[MemberReference, MemberInfo]):
    if isinstance(member, MemberInfo):
        print(f".. py:data:: {member.name}")
        if member.doc:
            print("\n" + prepend_spaces(member.doc.strip(), 3))
        print("\n")
    else:
        print(f".. py:data:: {member.name}")
        print(f"   :canonical: {member.qualified_name}")
        if member.doc:
            print("\n" + prepend_spaces(member.doc.strip(), 3))
        print("\n")


def annotation_to_str(annotation, prepend_tilde=True):
    if annotation and hasattr(annotation, "__name__"):
        name = full_class_name(annotation)
        if name.startswith("builtins."):
            name = name[9:]
            if name == "NoneType":
                return "None"
            return name
        if prepend_tilde:
            return "~" + name
        return name
    else:
        origin, args = get_origin(annotation), get_args(annotation)
        if origin is Union:
            return ' or '.join([annotation_to_str(arg) for arg in args])
        else:
            if prepend_tilde:
                return "~" + str(annotation)
            return str(annotation)


def generate_function_rst(member: Union[FunctionInfo, FunctionReference]):
    if isinstance(member, FunctionReference):
        canonical = f"{member.referee_module.qualified_name}.{member.name}"
    else:
        canonical = None

    print(f".. py:function:: {member.name}{format_signature(member.signature)}")
    with indent(3):
        if canonical:
            print(f":canonical: {canonical}")
        print("")
        if member.doc:
            print_docstring(member.doc)
            print("")

        sig = inspect.signature(member.method)
        for arg in sig.parameters.values():
            print(f":type {arg.name}: {annotation_to_str(arg.annotation)}")
        if sig.return_annotation:
            print(f":rtype: {annotation_to_str(sig.return_annotation)}")

    print("")


def print_docstring(doc):
    """
    Prints out a docstring and tries to correct formatting.

    This does two tasks:

    * Find the minimal indent of all the lines, and removes that from all lines.

    * Fixes indenting when using multiline definitions for :param ...: etc.

    * Inserts a space between :param: and non-blank lines
    """
    inset = None
    corrected_inset = None
    has_space = False
    for line in doc.splitlines():
        if len(line.strip()) > 0 and inset is None:
            inset = len(line) - len(line.lstrip(' '))
        if inset is None:
            continue
        if len(line.strip()) == 0:
            print("")
            corrected_inset = None
            has_space = True
        else:
            if line[inset:].startswith(":") and not has_space:
                print("")
            if corrected_inset and line[inset] != ":":
                print(" " * 4 + line.lstrip())
            else:
                print(line[inset:])
            has_space = False
        if line[inset:].startswith(":param") or line[inset:].startswith(":return"):
            corrected_inset = line[inset:].index(':', 1) + 2


def full_class_name(obj):
    if hasattr(obj, "__module__"):
        return f"{obj.__module__}.{obj.__name__}"
    return obj.__name__


def generate_tags(tag_dict: Dict[str, bool], html_class='sig-label'):
    tags = []
    for tag, valid in tag_dict.items():
        if valid:
            tags.append(f"<span class='{html_class}'>{tag}</span>")

    if tags:
        print(f".. raw:: html")
        with indent(3):
            print(f"")
            for tag in tags:
                print(f"{tag}")
            print(f"")


def generate_member_inheritance(member: Union[MemberInfo, PropertyInfo]):
    if member.inherits:
        for c, _ in member.inherits:
            print(
                f"Inherited from :py:obj:`{c.__name__} <{c.__module__}.{c.__name__}.{member.name}>`")
    if member.overrides:
        for c, _ in member.overrides:
            print(f"Overrides :py:obj:`{c.__name__} <{c.__module__}.{c.__name__}.{member.name}>`")
    if member.implements:
        for c, _ in member.implements:
            print(
                f"Implements :py:obj:`{c.__name__} <{c.__module__}.{c.__name__}.{member.name}>`")


def format_signature(sig: inspect.Signature):
    """Strip out annotations and remove self from a signature."""
    parameters = []
    for param_name, param_value in sig.parameters.items():
        if param_name == "self":
            continue
        new_param = inspect.Parameter(param_name, param_value.kind)
        parameters.append(new_param)
    return inspect.Signature(parameters)


def generate_class_method_rst(member: FunctionInfo):
    generate_tags({
        "Class Method": member.is_class_method,
        "Abstract": member.is_abstract,
        "Inherited": bool(member.inherits),
        "Overrides": bool(member.overrides),
        "Implementation": bool(member.implements)
    })

    if member.overloads:
        for overload in member.overloads:
            print(f".. py:method:: {member.name}{overload}")
            print("   :noindex:")
            print("   ")

    print(f".. py:method:: {member.name}{format_signature(member.signature)}")
    with indent(3):
        print("")
        if member.doc:
            print_docstring(member.doc)
            print("")


        sig = inspect.signature(member.method)
        for arg in sig.parameters.values():
            if arg.annotation != inspect._empty:
                print(f":type {arg.name}: {annotation_to_str(arg.annotation)}")
        if sig.return_annotation != inspect._empty:
            print(f":rtype: {annotation_to_str(sig.return_annotation)}")


        generate_member_inheritance(member)

    print("")

def generate_class_var_rst(member: MemberInfo):
    generate_tags({
        "Class Variable": member.is_classvar,
        "Inherited": bool(member.inherits),
        "Overrides": bool(member.overrides),
        "Implementation": bool(member.implements)
    })

    print(f"   .. py:attribute:: {member.name}")
    print("      ")
    if member.doc:
        print_docstring(member.doc, indent=6)
        print("      ")

    generate_member_inheritance(member)
    print("")


def generate_class_property_rst(prop: PropertyInfo):
    generate_tags({
        "Abstract": prop.is_abstract,
        "Event": prop.is_event,
        "Read Only": bool(prop.property.fget and not prop.property.fset),
        "Write Only": bool(prop.property.fset and not prop.property.fget),
        "Inherited": bool(prop.inherits),
        "Overrides": bool(prop.overrides),
        "Implementation": bool(prop.implements)
    })

    print(f".. py:property:: {prop.name}")
    with indent(3):
        if prop.is_abstract:
            print(f":abstractmethod:")
        sig = inspect.signature(prop.property.fget)
        if sig.return_annotation:
            print(f":type: {annotation_to_str(sig.return_annotation, prepend_tilde=False)}")
        print("")

        if prop.doc:
            print_docstring(prop.doc)
            print()

        generate_member_inheritance(prop)
        print("")

def generate_class_member_rst(member: Union[MemberInfo, FunctionInfo, PropertyInfo]):
    if isinstance(member, MemberInfo):
        generate_class_var_rst(member)
    elif isinstance(member, FunctionInfo):
        generate_class_method_rst(member)
    elif isinstance(member, PropertyInfo):
        generate_class_property_rst(member)
    else:
        raise ValueError(type(member))

def generate_class_rst(clss: ClassInfo):
    if isinstance(clss, ClassReference):
        clss = clss.target
        canonical = f"{clss.clss.__module__}.{clss.name}"
    else:
        canonical = None

    print(f".. py:class:: {clss.name}")

    with indent(3):

        if canonical:
            print(f":canonical: {canonical}")

        print("")
        if clss.doc:
            print_docstring(clss.doc)
            print("")

        if clss.protocols:
            print(
                f"Implements {', '.join([':py:class:`~' + full_class_name(protocol) + '`' for protocol in clss.protocols])}")
            print("")
        if clss.base_classes:
            print(
                f"Inherits from {', '.join([':py:class:`~' + full_class_name(base_class) + '`' for base_class in clss.base_classes])}")
            print("")

        class_magic_methods = {
            "with x": clss.has_method("__enter__") or clss.has_method("__exit__"),
            "len(x)": clss.has_method("__len__"),
            "abs(x)": clss.has_method("__abs__"),
            "round(x)": clss.has_method("__round__"),
            "~x": clss.has_method("__invert__"),
            "x + y": clss.has_method("__add__") or clss.has_method("__radd__"),
            "x - y": clss.has_method("__sub__") or clss.has_method("__rsub__"),
            "x * y": clss.has_method("__mul__") or clss.has_method("__rmul__"),
            "x @ y": clss.has_method("__matmul__") or clss.has_method("__rmatmul__"),
            "x << y": clss.has_method("__lshift__"),
            "x >> y": clss.has_method("__rshift__"),
            "x & y": clss.has_method("__and__"),
            "x | y": clss.has_method("__or__"),
            "x ^ y": clss.has_method("__xor__"),
            "x ** y": clss.has_method("__pow__"),
            "x[y]": clss.has_method("__getitem__"),
            "x[y] = z": clss.has_method("__setitem__"),
            "del x[y]": clss.has_method("__delitem__"),
            "y in x": clss.has_method("__contains__"),
            "for y in x": clss.has_method("__iter__"),
            "x(...)": clss.has_method("__call__"),
        }

        if any(class_magic_methods.values()):
            print("This class overrides the following operators:")
            print("")
            generate_tags(class_magic_methods, html_class="class-usage")
            print("")

        members = []

        for method in clss.methods:
            members.append(method)
        for property in clss.properties:
            members.append(property)
        for member in clss.members:
            members.append(member)

        for member in members:
            if member.name == "__init__":
                generate_class_member_rst(member)

        for member in members:
            if not member.name.startswith("_") and not member.inherits:
                generate_class_member_rst(member)


def generate_module_rst(module: ModuleInfo, folder: str):
    os.makedirs(folder, exist_ok=True)
    with open(f'{folder}/{module.qualified_name}.rst', 'w') as f, redirect_stdout(f):
        print(f"{module.qualified_name}")
        print("=" * len(module.qualified_name))
        print("\n")
        print(f".. py:module:: {module.qualified_name}")
        print("")
        print(f".. py:currentmodule:: {module.qualified_name}")
        print("")
        if module.doc:
            print("" + module.doc.strip())
            print("\n")

        for clss in module.classes:
            if not clss.name.startswith("_"):
                generate_class_rst(clss)

        for member in module.members:
            if not member.name.startswith("_"):
                generate_module_member_rst(member)

        for func in module.functions:
            if not func.name.startswith("_"):
                generate_function_rst(func)

        if module.modules:
            print(".. toctree::")
            print("   :maxdepth: 2")
            print("")

            for module in module.modules:
                if "_" in module.qualified_name:
                    continue
                generate_module_rst(module, folder)
                print(f"   {module.qualified_name}")
