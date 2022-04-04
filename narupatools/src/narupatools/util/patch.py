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

"""Techniques to patch existing functions."""

import types
from typing import Any


def compose(to_patch_func: Any) -> Any:
    """Alter an instance method of a class by compositing a second function."""
    if not isinstance(to_patch_func, types.MethodType):
        raise ValueError("Can't apply composition to non instance method.")

    def patcher(patching_func: Any) -> Any:
        def new_func(self: Any, *args: Any, **kwargs: Any) -> Any:
            return patching_func(self, to_patch_func(*args, **kwargs))

        dest_obj = to_patch_func.__self__
        dest_name = to_patch_func.__func__.__name__
        setattr(dest_obj, dest_name, types.MethodType(new_func, dest_obj))
        return patching_func

    return patcher


def replace(to_replace_func: Any) -> Any:
    """Replace the instance method of an object with a new function."""
    if not isinstance(to_replace_func, types.MethodType):
        raise ValueError("Can't apply composition to non instance method.")

    def replacer(replacing_func: Any) -> Any:
        def new_func(self: Any, *args: Any, **kwargs: Any) -> Any:
            return replacing_func(self, *args, **kwargs)

        dest_obj = to_replace_func.__self__
        dest_name = to_replace_func.__func__.__name__
        setattr(dest_obj, dest_name, types.MethodType(new_func, dest_obj))
        return replacing_func

    return replacer
