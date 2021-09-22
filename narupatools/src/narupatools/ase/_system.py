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
import numpy.typing as npt
from ase import Atoms
from infinite_sets import InfiniteSet
from narupa.trajectory import FrameData

from narupatools.core import UnitsNarupa
from narupatools.core.dynamics import DynamicsProperties, SimulationRotationProperties
from narupatools.frame._frame_source import FrameSource
from narupatools.physics._quaternion import quaternion
from narupatools.physics.typing import ScalarArray, Vector3Array, Vector3ArrayLike

from ._converter import ase_atoms_to_frame
from ._rotations import (
    get_angular_momenta,
    get_principal_moments,
    get_rotations,
    set_angular_momenta,
    set_angular_velocities,
    set_principal_moments,
    set_rotations,
)
from ._units import UnitsASE
from .calculators import NullCalculator
from ..override import override

_NarupaToASE = UnitsNarupa >> UnitsASE
_ASEToNarupa = UnitsASE >> UnitsNarupa


def create_ase_atoms(*, positions: Vector3ArrayLike, velocities: Vector3ArrayLike, masses: ScalarArray):
    """
    Create a new ASE atoms object in Narupa units.

    :param positions: Position of each atom in nanometers.
    :param velocities: Velocity of each atom in nanometers per picoseconds.
    :param masses: Mass of each atom in daltons.
    :return: ASE atoms object.
    """
    kwargs = {}

    if positions is not None:
        kwargs["positions"] = positions * _NarupaToASE.length

    if velocities is not None:
        kwargs["velocities"] = velocities * _NarupaToASE.velocity

    if masses is not None:
        kwargs["masses"] = masses * _NarupaToASE.mass

    atoms = Atoms(**kwargs)

    atoms.calc = NullCalculator()

    return atoms


class ASESystem(FrameSource, DynamicsProperties, SimulationRotationProperties):
    """Wrapper around an ASE `Atoms` object so it is exposed consistently."""

    def __init__(self, atoms: Atoms):
        """
        Create a wrapper around the given ASE `Atoms` object.

        :param atoms: ASE `Atoms` to wrap.
        """
        self._atoms = atoms

    @override
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
    @override
    def positions(self) -> Vector3Array:  # noqa: D102
        return self.atoms.positions * _ASEToNarupa.length  # type: ignore

    @positions.setter
    def positions(self, value: Vector3ArrayLike) -> None:
        self.atoms.set_positions(np.asfarray(value) * _NarupaToASE.length)

    @property
    @override
    def velocities(self) -> Vector3Array:  # noqa: D102
        return self.atoms.get_velocities() * _ASEToNarupa.velocity  # type: ignore

    @velocities.setter
    def velocities(self, value: Vector3ArrayLike) -> None:
        self.atoms.set_velocities(np.asfarray(value) * _NarupaToASE.velocity)

    @property
    @override
    def forces(self) -> Vector3Array:  # noqa: D102
        return self.atoms.get_forces() * _ASEToNarupa.force  # type: ignore

    @property
    @override
    def masses(self) -> ScalarArray:  # noqa: D102
        return self.atoms.get_masses() * _ASEToNarupa.mass  # type: ignore

    @masses.setter
    def masses(self, value: ScalarArray) -> None:
        self.atoms.set_masses(np.asfarray(value) * _NarupaToASE.mass)

    @property
    @override
    def kinetic_energy(self) -> float:  # noqa: D102
        return self.atoms.get_kinetic_energy() * _ASEToNarupa.energy

    @override
    @property
    def potential_energy(self) -> float:  # noqa: D102
        return self.atoms.get_potential_energy() * _ASEToNarupa.energy

    @property
    @override
    def orientations(self) -> npt.NDArray[quaternion]:  # noqa: D102
        return get_rotations(self.atoms)

    @orientations.setter
    def orientations(self, value: npt.NDArray[quaternion]) -> None:
        set_rotations(self.atoms, value)

    @property
    @override
    def angular_momenta(self) -> Vector3Array:  # noqa: D102
        return get_angular_momenta(self.atoms)

    @angular_momenta.setter
    def angular_momenta(self, value: Vector3ArrayLike) -> None:
        set_angular_momenta(
            self.atoms, np.asfarray(value) * _ASEToNarupa.angular_momenta
        )

    @property
    @override
    def angular_velocities(self) -> Vector3Array:
        """
        Angular velocity of each particle abouts its center of mass.

        :return: Array of angular velocities in radians per picoseconds.
        """
        raise AttributeError

    @angular_velocities.setter
    def angular_velocities(self, value: Vector3ArrayLike) -> None:
        set_angular_velocities(
            self.atoms, np.asfarray(value) * _NarupaToASE.angular_velocity
        )

    @property
    @override
    def moments_of_inertia(self) -> Vector3Array:  # noqa: D102
        return get_principal_moments(self.atoms) * _ASEToNarupa.moment_inertia

    @moments_of_inertia.setter
    def moments_of_inertia(self, value: Vector3ArrayLike) -> None:
        set_principal_moments(
            self.atoms, np.asfarray(value) * _NarupaToASE.moment_inertia
        )