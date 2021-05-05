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

from typing import Protocol

import pytest

from narupatools.core.event import CallbackMissingParametersError, Event


class CallbackWithOneArgument(Protocol):
    def __call__(self, *, arg1):
        pass


def test_event_add_callback_raises_warning():
    def callback_without_kwargs(*, arg1):
        pass

    event = Event()

    with pytest.warns(UserWarning):
        event.add_callback(callback_without_kwargs)


def test_event_add_callback_with_wrong_signature():
    def callback_right_signature(*, arg1):
        pass

    def callback_wrong_signature(*, arg2):
        pass

    def callback_wrong_signature_with_kwargs(*, arg2, **kwargs):
        pass

    def callback_kwargs(**kwargs):
        pass

    event = Event(CallbackWithOneArgument)

    with pytest.raises(CallbackMissingParametersError):
        event.add_callback(callback_wrong_signature)

    event.add_callback(callback_right_signature)

    event.add_callback(callback_wrong_signature_with_kwargs)

    event.add_callback(callback_kwargs)
