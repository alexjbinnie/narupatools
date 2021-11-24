import functools
import inspect
import sys
import warnings
from contextlib import contextmanager
from typing import Protocol

from narupatools.override import marked_as_override


def get_src_link(class_, name_, item_):
    try:
        _, lineno = inspect.getsourcelines(item_)
        src = f"{inspect.getsourcefile(item_)}:{lineno} {name_}"
    except (ValueError, TypeError):
        _, lineno = inspect.getsourcelines(class_)
        src = f"{inspect.getsourcefile(class_)}:{lineno} {name_}"
    return src


def isabstract(obj):
    return hasattr(obj, "__isabstractmethod__")


def accepts_variable_positionals(sig):
    return any(param.kind == inspect.Parameter.VAR_POSITIONAL for param in sig.parameters.values())


def accepts_variable_keywords(sig):
    return any(param.kind == inspect.Parameter.VAR_KEYWORD for param in sig.parameters.values())


def compare_signatures(item, original):
    """Check if the first method is a valid override of the second method."""
    sig_ovr = inspect.signature(item)
    try:
        if isinstance(original, property):
            sig_orig = inspect.signature(original.fget)
        else:
            sig_orig = inspect.signature(original)
    except ValueError:
        return

    params_orig = list(sig_orig.parameters.items())
    params_ovr = list(sig_ovr.parameters.items())

    for index, params in enumerate(sig_orig.parameters.items()):
        name, param = params
        if param.kind == inspect.Parameter.POSITIONAL_OR_KEYWORD:
            try:
                name2, param2 = params_ovr[index]
                if name != name2:
                    warnings.warn(f"Parameter name {name2} of {item} incompatible with {name} of {original}", UserWarning)
                if param2.kind == inspect.Parameter.POSITIONAL_ONLY:
                    warnings.warn(
                        f"Overriden function {item} has parameter {name2} as positional, whilst {original} allows it as a keyword.",
                        UserWarning)
            except IndexError:
                warnings.warn(f"Parameter name {name} of {item} missing in override.", UserWarning)

        if param.kind == inspect.Parameter.POSITIONAL_ONLY:
            _, param2 = params_ovr[index]
            if param2.kind not in (inspect.Parameter.POSITIONAL_ONLY, inspect.Parameter.POSITIONAL_OR_KEYWORD):
                warnings.warn(
                    f"Overriden function {item} does not support positional-only argument {name} from {original}",
                    UserWarning)

        if param.kind == inspect.Parameter.KEYWORD_ONLY:
            if param.name not in sig_ovr.parameters and not accepts_variable_keywords(sig_ovr):
                warnings.warn(f"Overriden function {item} does not support keyword {param.name} from {original}",
                              UserWarning)
            if name in sig_ovr.parameters and sig_ovr.parameters[name].kind not in (
                    inspect.Parameter.POSITIONAL_OR_KEYWORD, inspect.Parameter.KEYWORD_ONLY):
                warnings.warn(
                    f"Overriden function {item} does not support {name} as a keyword argument, unlike {original}.",
                    UserWarning)

    for index, params in enumerate(params_ovr):
        name, param = params
        if param.kind == inspect.Parameter.POSITIONAL_ONLY:
            if len(params_orig) < index + 1:
                warnings.warn(f"{item} has more parameters than {original}",
                              UserWarning)
            else:
                _, param_orig = params_orig[index]
                if param_orig.kind != inspect.Parameter.POSITIONAL_ONLY:
                    warnings.warn(
                        f"{item} has positional-only parameter {name} which is not positional only in {original}",
                        UserWarning)

    if accepts_variable_positionals(sig_orig) and not accepts_variable_positionals(sig_ovr):
        warnings.warn(
            f"Overriden function {item} does not accept a variable number of positional arguments, but the original function {original} does.",
            UserWarning)
    if accepts_variable_keywords(sig_orig) and not accepts_variable_keywords(sig_ovr):
        warnings.warn(
            f"Overriden function {item} does not accept a variable number of keyword arguments, but the original function {original} does.",
            UserWarning)


def explore_inheritance(class_, name, item, info):
    for cls in inspect.getmro(class_):
        if cls is class_:
            continue
        if hasattr(cls, name):
            original = getattr(cls, name)
            if original is item or (
                    inspect.ismethod(item) and inspect.ismethod(original) and item.__func__ is getattr(cls,
                                                                                                       name).__func__):
                info.inherits.append((cls, original))
                return
            if not marked_as_override(item) and (
                    inspect.ismethod(item) or inspect.isfunction(item) or isinstance(item,
                                                                                     property)) and not name.startswith(
                "__"):
                warnings.warn(f"{get_src_link(class_, name, item)} needs @override", UserWarning)
            overriding_abstract = False
            if inspect.ismethod(item) or inspect.isfunction(item):
                if not name.startswith("__"):
                    compare_signatures(item, original)
                if isabstract(original):
                    overriding_abstract = True
            if overriding_abstract:
                info.implements.append((cls, original))
            else:
                info.overrides.append((cls, original))
            return
    if marked_as_override(item):
        warnings.warn(
            f"{get_src_link(class_, name, item)} marked with @override but does not override anything. {class_.__name__} {name}",
            UserWarning)


def is_class_method(class_, name_, func_):
    return inspect.ismethod(func_) and func_.__self__ is class_


def prepend_tab(s):
    lines = s.splitlines(keepends=True)
    return "".join(["\t" + line for line in lines])


def prepend_spaces(s, n=3):
    lines = s.splitlines(keepends=True)
    return "".join([(" " * n) + line for line in lines])


@contextmanager
def indent(indent=4, space=chr(32)):
    """ Precede each carriage return with some quantity of spaces.
        While nesting `indented_output` contexts, be prepared
        that if a line does not begin with a line feed character,
        it shall print with unchanged indentation.
        To keep the indentation as expected,
        print last line feed character of a parent block
        in the beginning of a child block.
    """
    try:
        original_stderr_write = sys.stderr.write
        original_stdout_write = sys.stdout.write

        def write_indented(write_call, data):
            write_call.__call__((space * indent) + data)

        sys.stderr.write = functools.partial(
            write_indented, original_stderr_write)
        sys.stdout.write = functools.partial(
            write_indented, original_stdout_write)

        yield
    finally:
        sys.stderr.write = original_stderr_write
        sys.stdout.write = original_stdout_write


import collections.abc as collections_abc

collections_abc = [getattr(collections_abc, name) for name in collections_abc.__all__]


def is_protocol(clss):
    return issubclass(clss, Protocol) or clss in collections_abc
