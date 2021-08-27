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

import numpy as np
from ase import Atoms
from infinite_sets import InfiniteSet
from narupa.trajectory import FrameData

from narupatools.core import UnitsNarupa
from narupatools.core.dynamics import DynamicsProperties
from narupatools.frame._frame_source import FrameSource
from narupatools.physics.typing import ScalarArray, Vector3Array, Vector3ArrayLike

from ._converter import ase_atoms_to_frame
from ._units import UnitsASE
from .calculators import NullCalculator

_NarupaToASE = UnitsNarupa >> UnitsASE
_ASEToNarupa = UnitsASE >> UnitsNarupa


class ASESystem(FrameSource, DynamicsProperties):
    """Wrapper around an ASE `Atoms` object so it is exposed consistently."""

    def __init__(self, atoms: Atoms):
        """
        Create a wrapper around the given ASE `Atoms` object.

        :param atoms: ASE `Atoms` to wrap.
        """
        self._atoms = atoms

    def get_frame(self, fields: InfiniteSet[str]) -> FrameData:  # noqa: D102
        frame = FrameData()
        ase_atoms_to_frame(self._atoms, fields=fields, frame=frame)
        return frame

    @property
    def atoms(self) -> Atoms:
        """Wrapped ASE atoms object represented by this system."""
        return self._atoms

    @staticmethod
    def create(atoms: Atoms) -> ASESystem:
        """
        Create a an ASE system which can be broadcast through a session.

        If the Atoms object does not have a calculator, a null calculator will be added.

        :param atoms: ASE atoms object.
        :return: Representation of the system which can be broadcast through a session.
        """
        if atoms.get_calculator() is None:
            atoms.set_calculator(NullCalculator())
        return ASESystem(atoms)

    @property
    def positions(self) -> Vector3Array:  # noqa: D102
        return self.atoms.positions * _ASEToNarupa.length  # type: ignore

    @positions.setter
    def positions(self, value: Vector3ArrayLike) -> None:
        self.atoms.set_positions(np.asfarray(value) * _NarupaToASE.length)

    @property
    def velocities(self) -> Vector3Array:  # noqa: D102
        return self.atoms.get_velocities() * _ASEToNarupa.velocity  # type: ignore

    @velocities.setter
    def velocities(self, value: Vector3ArrayLike) -> None:
        self.atoms.set_velocities(np.asfarray(value) * _NarupaToASE.velocity)

    @property
    def forces(self) -> Vector3Array:  # noqa: D102
        return self.atoms.get_forces() * _ASEToNarupa.force  # type: ignore

    @property
    def masses(self) -> ScalarArray:  # noqa: D102
        return self.atoms.get_masses() * _ASEToNarupa.mass  # type: ignore

    @property
    def kinetic_energy(self) -> float:  # noqa: D102
        return self.atoms.get_kinetic_energy() * _ASEToNarupa.energy

    @property
    def potential_energy(self) -> float:  # noqa: D102
        return self.atoms.get_potential_energy() * _ASEToNarupa.energy
