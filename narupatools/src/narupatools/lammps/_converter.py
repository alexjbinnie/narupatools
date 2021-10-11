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

"""Converter functions for interfacing with LAMMPS."""
from typing import Any, Optional, Type, TypeVar, Union

from ase.atoms import Atoms
from infinite_sets import InfiniteSet
from narupa.trajectory import FrameData

from narupatools.frame import FrameConverter, convert

from ._calculator import LAMMPSCalculator
from ._simulation import LAMMPSSimulation
from ..override import override

_TType = TypeVar("_TType")


class LAMMPSConverter(FrameConverter):
    """FrameConverter for the mdtraj package."""

    @classmethod
    @override(FrameConverter.convert_to_frame)
    def convert_to_frame(  # noqa: D102
        cls, object_: Any, /, *, fields: InfiniteSet[str], existing: Optional[FrameData]
    ) -> FrameData:
        if isinstance(object_, LAMMPSSimulation):
            return object_.get_frame(fields=fields, existing=existing)
        raise NotImplementedError

    @classmethod
    @override(FrameConverter.convert_from_frame)
    def convert_from_frame(  # noqa: D102
        cls,
        frame: FrameData,
        destination: Union[Type[_TType], _TType],
        *,
        fields: InfiniteSet[str],
    ) -> _TType:
        raise NotImplementedError


def atoms_from_lammps_simulation(simulation: LAMMPSSimulation) -> Atoms:
    """
    Create an ASE Atoms object based on a LAMMPS simulation.

    :param simulation: LAMMPS simulation to use.
    :return: ASE atoms object with a calculator that references the given simulation.
    """
    atoms = convert(simulation, Atoms)
    calc = LAMMPSCalculator(simulation, atoms)
    atoms.set_calculator(calc)
    return atoms
