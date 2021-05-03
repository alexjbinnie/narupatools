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

"""Simple wrapper around an ASE atoms object."""

from __future__ import annotations

from typing import Optional

from MDAnalysis import Universe
from ase import Atoms
from infinite_sets import InfiniteSet
from narupa.trajectory import FrameData

from narupatools.ase import NullCalculator
from narupatools.ase import ase_atoms_to_frame
from narupatools.frame.frame_source import FrameSource
from narupatools.mdanalysis import mdanalysis_universe_to_frame


class ASESystem(FrameSource):
    """Wrapper around an ASE `Atoms` object so it is exposed consistently."""

    def __init__(self, atoms: Atoms, *, universe: Optional[Universe] = None):
        """
        Create a wrapper around the given ASE `Atoms` object.

        :param atoms: ASE `Atoms` to wrap.
        :param universe: Optional MDAnalysis universe containing additional information.
        """
        self._atoms = atoms
        self._universe = universe

    def get_frame(self, fields: InfiniteSet[str]) -> FrameData:  # noqa: D102
        frame = FrameData()
        if self._universe:
            ase_atoms_to_frame(self._atoms, fields=fields, frame=frame)
            added_fields = set(frame.arrays.keys()) | set(frame.values.keys())
            mdanalysis_universe_to_frame(
                self._universe, fields=fields - added_fields, frame=frame
            )
            return frame
        else:
            ase_atoms_to_frame(self._atoms, fields=fields, frame=frame)
        return frame

    @property
    def atoms(self) -> Atoms:
        """Wrapped ASE atoms object represented by this system."""
        return self._atoms

    @staticmethod
    def create(atoms: Atoms, *, universe: Optional[Universe] = None) -> ASESystem:
        """
        Create a an ASE system which can be broadcast through a session.

        If the Atoms object does not have a calculator, a null calculator will be added.

        :param atoms: ASE atoms object.
        :param universe: Optional MDAnalysis universe with additional system
                         information.
        :return: Representation of the system which can be broadcast through a session.
        """
        if atoms.get_calculator() is None:
            atoms.set_calculator(NullCalculator())
        return ASESystem(atoms, universe=universe)
