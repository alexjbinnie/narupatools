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

"""Provides simple decorator for marking methods as overrides."""

from typing import Any, TypeVar

_T = TypeVar("_T")


def override(f: _T, /) -> _T:
    """
    Mark a method or property as overriding a method in a base class.

    This is a decorator that can be applied to methods and properties:

    .. code-block:: python

       class MyBaseClass:

           def my_method():
               ...

       class MyClass(MyBaseClass):

           @override
           def my_method():
               ...

    This does not perform any checks to ensure that there is actually a method to override. This
    merely annotates the method or property.

    To check if a method is overriden, use :obj:`marked_as_override`.
    """
    if isinstance(f, property):
        if f.fget is not None:
            f.fget.__override__ = True  # type: ignore[attr-defined]
        if f.fset is not None:
            f.fset.__override__ = True  # type: ignore[attr-defined]
        if f.fdel is not None:
            f.fdel.__override__ = True  # type: ignore[attr-defined]
    else:
        f.__override__ = True  # type: ignore[attr-defined]
    return f


def marked_as_override(f: Any) -> bool:
    """
    Check if the given method or property is marked using the @override annotation.

    This merely checks if the argument has been annotated with @override. It does
    not actually check if the method or property overrides something in the base class.
    """
    if isinstance(f, property):
        return hasattr(f.fget, "__override__")
    else:
        return hasattr(f, "__override__")
