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

"""show_*** commands for showing an NGLWidget for a given system."""

from typing import Any

from ase.atoms import Atoms
from narupa.trajectory import FrameData
from nglview import NGLWidget

from ..frame.frame_source import TrajectorySource
from .structure import ASEStructure, FrameDataStructure, NarupaToolsTrajectory


def show_ase(atoms: Atoms, /, **kwargs: Any) -> NGLWidget:
    """
    Open a NGLWidget showing the given ASE atoms object.

    This works in the same was as `nglview.show_ase`, except that
    instead of writing to a temporary file it writes to a string in memory. This
    prevents permission errors if the temporary file cannot be written (for example, in
    Jupyter notebook).

    :param atoms: An ASE Atoms object to show.
    :param kwargs: Arguments for NGLWidget.
    :return: An NGLWidget showing the provided atoms.
    """
    structure = ASEStructure(atoms)
    return NGLWidget(structure, **kwargs)


def show_narupa(frame_data: FrameData, /, **kwargs: Any) -> NGLWidget:
    """
    Open a NGLWidget showing the given Narupa `FrameData`.

    :param frame_data: A Narupa `FrameData` to visualize.
    :param kwargs: Arguments for NGLWidget.
    :return: An NGLWidget showing the provided frame data.
    """
    structure = FrameDataStructure(frame_data)
    return NGLWidget(structure, **kwargs)


def show_narupatools_traj(source: TrajectorySource, /, **kwargs: Any) -> NGLWidget:
    """
    Open a NGLWidget showing the given narupatools `TrajectorySource`.

    :param source: Source to visualize.
    :param kwargs: Arguments for NGLWidget.
    :return: An NGLWidget showing the provided object.
    """
    structure = NarupaToolsTrajectory(source)
    return NGLWidget(structure, **kwargs)
