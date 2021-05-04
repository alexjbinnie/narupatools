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
import os
from contextlib import contextmanager

from narupatools.core.event import EventListener


@contextmanager
def assert_event_called(event: EventListener):
    was_called = False

    def on_called(*args, **kwargs):
        nonlocal was_called
        was_called = True

    event.add_callback(on_called)

    yield

    event.remove_callback(on_called)

    assert was_called


@contextmanager
def assert_event_not_called(event: EventListener):
    was_called = False

    def on_called(*args, **kwargs):
        nonlocal was_called
        was_called = True

    event.add_callback(on_called)

    yield

    event.remove_callback(on_called)

    assert not was_called


def add_mark(*, filename, mark, items):
    dir_path = os.path.dirname(os.path.realpath(filename))
    for item in items:
        file_path = os.path.realpath(item.fspath)
        if file_path.startswith(dir_path):
            item.add_marker(mark)
