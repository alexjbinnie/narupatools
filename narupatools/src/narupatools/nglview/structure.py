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

"""Classes derived from Structure that allow NGLView to render new systems."""

from io import StringIO

import nglview
from ase.atoms import Atoms
from ase.io import write
from narupa.trajectory.frame_data import FrameData

from narupatools.frame.converter import frame_to_pdb_string


class ASEStructure(nglview.Structure):
    """ASE Structure for nglview that does not use temporary files."""

    def __init__(self, atoms: Atoms, /) -> None:
        super().__init__()
        self.atoms = atoms

    def get_structure_string(self) -> str:
        """Create a PDB string so NGLView can read in the structure."""
        file = StringIO("")
        write(file, self.atoms, "proteindatabank")
        return file.getvalue()


class FrameDataStructure(nglview.Structure):
    """FrameData Structure for nglview that does not use temporary files."""

    def __init__(self, frame: FrameData, /) -> None:
        super().__init__()
        self.frame = frame

    def get_structure_string(self) -> str:
        """Create a PDB string so NGLView can read in the structure."""
        return frame_to_pdb_string(self.frame)
