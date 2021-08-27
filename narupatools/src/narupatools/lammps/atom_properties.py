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

"""Hard-coded atom properties defined by LAMMPS."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Generic, Literal, TypeVar, Union, overload

import numpy as np
import numpy.typing as npt

from ._constants import VariableType
from .exceptions import LAMMPSError

_TReturnType = TypeVar("_TReturnType")


@dataclass
class AtomProperty(Generic[_TReturnType]):
    """
    Definition of an atom property defined by LAMMPS.

    These definitions can be found in `atom.cpp` in the LAMMPS src/ folder.
    """

    key: str
    """Key used to identify property in LAMMPS, such as 'x' or 'f'."""
    type: VariableType
    """Enum which informs LAMMPS what data type the property is internally."""
    components: int
    """Number of components for each atom value."""

    @overload
    @classmethod
    def define(
        cls,
        key: str,
        type: Literal[VariableType.INTEGER, VariableType.INTEGER_ARRAY],
        components: int,
    ) -> AtomProperty[npt.NDArray[np.int64]]:
        ...

    @overload
    @classmethod
    def define(
        cls,
        key: str,
        type: Literal[VariableType.DOUBLE, VariableType.DOUBLE_ARRAY],
        components: int,
    ) -> AtomProperty[npt.NDArray[np.float64]]:
        ...

    @classmethod
    def define(
        cls,
        key: str,
        type: VariableType,
        components: int,
    ) -> Union[
        AtomProperty[npt.NDArray[np.float64]], AtomProperty[npt.NDArray[np.int64]]
    ]:
        """
        Define an atom property used by LAMMPS.

        This property can be used with the `gather_atoms` function to collect the property for all atoms,
        ordered by their atom ID.

        Calling this should only be necessary if a new custom property has been implemented in LAMMPS, as
        all the default properties are already defined in :mod:`narupatools.lammps.atom_properties`

        :param key: Key used to identify property in LAMMPS, such as 'x' or 'f'.
        :param type: The type of the returned value, either a floating point number or an integer.
        :param components: The number of components in this property. Normally '1' for scalar values such as
                          charge, and '3' for vector values such as position.
        :return: An :class:`AtomProperty` that can be used with :meth:`LAMMPSSimulation.gather_atoms`.
        """
        if type not in [
            VariableType.INTEGER,
            VariableType.INTEGER_ARRAY,
            VariableType.DOUBLE,
            VariableType.DOUBLE_ARRAY,
        ]:
            raise LAMMPSError(f"Atom property cannot be of variable type {type}")
        return AtomProperty(key, type, components)  # type: ignore[return-value]


AtomID = AtomProperty.define("id", VariableType.INTEGER, 1)

AtomType = AtomProperty.define("type", VariableType.INTEGER, 1)
"""Integer atom type."""

GroupMask = AtomProperty.define("mask", VariableType.INTEGER, 1)
"""Bitmask defining membership of each LAMMPS group."""

ImageFlags = AtomProperty.define("image", VariableType.INTEGER, 1)
"""Encoded flag defining what image each atom is in."""

Position = AtomProperty.define("x", VariableType.DOUBLE_ARRAY, 3)
"""Position of each atom."""

Velocity = AtomProperty.define("v", VariableType.DOUBLE_ARRAY, 3)
"""Velocity of each atom."""

Force = AtomProperty.define("f", VariableType.DOUBLE_ARRAY, 3)
"""Force on each atom."""

MoleculeID = AtomProperty.define("molecule", VariableType.INTEGER, 1)
"""ID of the molecule each atom belongs to."""

Charge = AtomProperty.define("q", VariableType.DOUBLE, 1)
"""Charge of each atom."""

DipoleMoment = AtomProperty.define("mu", VariableType.DOUBLE_ARRAY, 3)
"""Dipole moment of each atom."""

AngularVelocity = AtomProperty.define("omega", VariableType.DOUBLE_ARRAY, 3)
"""Angular velocity of each atom."""

AngularMomentum = AtomProperty.define("angmom", VariableType.DOUBLE_ARRAY, 3)
"""Angular momentum of each atom."""

Torque = AtomProperty.define("torque", VariableType.DOUBLE_ARRAY, 3)
"""Torque on each atom."""

Radius = AtomProperty.define("radius", VariableType.DOUBLE, 1)
"""Per atom radius."""

PerAtomMass = AtomProperty.define("rmass", VariableType.DOUBLE, 1)
"""Per atom mass."""

EllipsoidParticle = AtomProperty.define("ellipsoid", VariableType.INTEGER, 1)
"""Flag whether each atom is an Ellipsoid or not."""

LineParticle = AtomProperty.define("line", VariableType.INTEGER, 1)
"""Flag whether each atom is a Line particle or not."""

TriangulatedParticle = AtomProperty.define("tri", VariableType.INTEGER, 1)
"""Flag whether each atom is a Triangulated particle or not."""

BodyParticle = AtomProperty.define("body", VariableType.INTEGER, 1)
"""Flag whether each atom is a Body particle or not."""
