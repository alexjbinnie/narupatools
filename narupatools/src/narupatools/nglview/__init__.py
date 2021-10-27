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

"""Code for interfacing with NGLView."""

import importlib
from typing import Any

from ase import Atoms
from narupa.trajectory import FrameData
from nglview import NGLWidget

from ..app import Client, Session
from ..core.dynamics import SimulationDynamics

__has_ngl = importlib.util.find_spec("nglview") is not None

if not __has_ngl:
    raise ImportError("narupatools.nglview requires nglview to be installed.")

from ._client import show_client
from ._dynamics import show_dynamics
from ._session import show_session
from ._show import show_ase, show_narupa, show_trajectory
from ._structure import ASEStructure, FrameDataStructure, FrameDataTrajectory


def show(obj: Any) -> NGLWidget:
    if isinstance(obj, Atoms):
        return show_ase(obj)
    if isinstance(obj, Client):
        return show_client(obj)
    if isinstance(obj, Session):
        return show_session(obj)
    if isinstance(obj, SimulationDynamics):
        return show_dynamics(obj)
    if isinstance(obj, FrameData):
        return show_narupa(obj)
    if isinstance(obj, list) and isinstance(obj[0], FrameData):
        return NGLWidget(FrameDataTrajectory(obj))
    raise ValueError(f"Cannot work out how to show object {obj} using nglview")


__all__ = [
    "show_ase",
    "show_client",
    "show_narupa",
    "show_dynamics",
    "show_trajectory",
    "show_session",
    "ASEStructure",
    "FrameDataStructure",
]
