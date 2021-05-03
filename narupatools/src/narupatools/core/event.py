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
#
# Originally part of the narupa-core package.
# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Modified under the terms of the GPL.

"""Event-based callback system with support for typing."""

from __future__ import annotations

import contextlib
import warnings
from inspect import Parameter, Signature, signature
from typing import (
    AbstractSet,
    Any,
    Callable,
    Generic,
    Iterable,
    List,
    Optional,
    Protocol,
    Type,
    TypeVar,
)

TCallback = TypeVar("TCallback", contravariant=True)


class CallbackMissingParametersError(SyntaxError):
    """Error thrown when a callback is missing a named parameter."""

    def __init__(self, missing_params: Iterable[str], /):
        super().__init__(
            f"Callback missing parameter(s) with name(s) {missing_params}. "
            f"Either add these parameters to the callback or use **kwargs."
        )


class EventListener(Protocol[TCallback]):
    """Protocol describing an event as seen from an external class."""

    def add_callback(self, callback: TCallback) -> None:
        """
        Add a callback.

        :param callback: Function to be invoked when the event is triggered.
        """
        ...

    def remove_callback(self, callback: TCallback) -> None:
        """
        Remove a callback.

        :param callback: Function previous added using
                         :ref:`EventListener.add_callback`.
        """
        ...


TEventCallback = TypeVar("TEventCallback", bound=Callable[..., None])


class Event(EventListener, Generic[TEventCallback]):
    """
    Class which stores a set of callbacks, which are invoked when an event is published.

    The narupatools version of Event performs some simple type checking, to ensure that
    callbacks support all the parameters that this event produces. It also prints a
    warning if an added callback does not handle generic ``**kwargs``, which is good
    practise that future proof callbacks.
    """

    _callbacks: List[TEventCallback]
    _parameter_names: AbstractSet[str]

    def __init__(self, callback_type: Optional[Type[TEventCallback]] = None) -> None:
        r"""
        Create a new event, optionally passing a type that represents the callback.

        :param callback_type: Desired type of the callback. If this has a ``__call__``
                              method, the named arguments of this will be stored. Adding
                              a callback to this event which cannot handle all these
                              named parameters (either directly or by using
                              ``**kwargs``) will raise an exception.
        """
        self._callbacks = []
        self._parameter_names = set()
        if callback_type is not None:
            with contextlib.suppress(AttributeError):
                self._parameter_names = _get_named_parameters(
                    signature(callback_type.__call__)
                )

    def add_callback(self, callback: TEventCallback) -> None:
        """
        Add a callback to this event.

        :param callback: The callback to be called when this event is triggered.
        :raises CallbackMissingParametersError: Callback is missing one or more keyword
                                                arguments the event uses.
        """
        callback_signature = signature(callback)

        has_kwargs = any(
            param.kind == Parameter.VAR_KEYWORD
            for param in callback_signature.parameters.values()
        )

        if not has_kwargs:
            warnings.warn(
                "Adding callback which doesn't take **kwargs. "
                "This makes it susceptible to breaking if new parameters are added to "
                "the event."
            )

        if not has_kwargs and len(self._parameter_names) > 0:
            params = _get_named_parameters(callback_signature)
            missing_params = self._parameter_names - params
            if len(missing_params) > 0:
                raise CallbackMissingParametersError(missing_params)

        self._callbacks.append(callback)

    def remove_callback(self, callback: TEventCallback) -> None:
        """
        Remove a callback from this event.

        This operation is atomic, and will not fail if ``callback`` is not listening to
        this event.

        :param callback: The callback to be removed from this event's callbacks
        """
        self._callbacks.remove(callback)

    @property
    def invoke(self) -> TEventCallback:
        """
        Invoke the callbacks associated with this event with the provided arguments.

        This is implemented as a property which returns the function that should then be
        invoked. This is so the invoke function can be correctly typed, and allows IDE's
        to autocomplete arguments correctly.
        """

        def _invoke(*args: Any, **kwargs: Any) -> None:
            for callback in self._callbacks:
                callback(*args, **kwargs)

        return _invoke  # type: ignore


def _get_named_parameters(func_signature: Signature) -> AbstractSet[str]:
    return {
        name
        for name, param in func_signature.parameters.items()
        if (
            param.kind == Parameter.POSITIONAL_OR_KEYWORD
            or param.kind == Parameter.KEYWORD_ONLY
        )
        and name != "self"
    }
