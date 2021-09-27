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

"""Regions in LAMMPS simulations."""

from __future__ import annotations

from abc import ABCMeta, abstractmethod
from typing import Any

from narupatools.physics.typing import Vector3
from narupatools.physics.units import UnitConversion
from narupatools.physics.vector import vector


class Region:
    """Region defined in a LAMMPS simulation."""

    def __init__(
        self, simulation: Any, region_id: str, specification: RegionSpecification
    ):
        self._simulation = simulation
        self._id = region_id
        self._specification = specification

    @property
    def region_id(self) -> str:
        """ID of this region."""
        return self._id

    @property
    def specification(self) -> RegionSpecification:
        """Specification used for this region."""
        return self._specification


class RegionSpecification(metaclass=ABCMeta):
    """Specification for a type of LAMMPS region."""

    @property
    @abstractmethod
    def style(self) -> str:
        """Style of the region."""

    @abstractmethod
    def args(self, conversion: UnitConversion) -> str:
        """
        Parameters for the region.

        :param conversion: Unit conversion from Narupa to the units of the given
                           simulation.
        """


class Box(RegionSpecification):
    """Specification of an axis-aligned box for a LAMMPS simulation."""

    def __init__(
        self, xlo: float, xhi: float, ylo: float, yhi: float, zlo: float, zhi: float
    ):
        self.lo = vector(xlo, ylo, zlo)
        self.hi = vector(xhi, yhi, zhi)

    @classmethod
    def bounds(cls, lower: Vector3, upper: Vector3) -> Box:
        """
        Create a specification for an axis-aligned box for a LAMMPS simulation.

        :param lower: Lower corner (xlo, ylo, zlo).
        :param upper: Upper corner (xhi, yhi, zhi).
        :return: Box specification with given bounds.
        """
        return Box(
            xlo=lower[0],
            xhi=upper[0],
            ylo=lower[1],
            yhi=upper[1],
            zlo=lower[2],
            zhi=upper[2],
        )

    @classmethod
    def from_size(cls, *, size: Vector3, center: Vector3 = vector(0, 0, 0)) -> Box:
        """
        Create a specification for an axis-aligned box for a LAMMPS simulation.

        :param size: Box sides in nanometers..
        :param center: Box center in nanometers..
        :return: Box specification with given size and center.
        """
        return cls.bounds(lower=center - size / 2.0, upper=center + size / 2.0)

    @property
    def style(self) -> str:  # noqa: D102
        return "block"

    def args(self, conversion: UnitConversion) -> str:  # noqa: D102
        lo = self.lo * conversion.length
        hi = self.hi * conversion.length
        return f"{lo[0]} {hi[0]} {lo[1]} {hi[1]} {lo[2]} {hi[2]}"
