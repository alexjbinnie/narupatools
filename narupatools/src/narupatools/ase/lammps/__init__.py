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

"""Integration between LAMMPS and ASE."""

from ase import Atoms

from narupatools.frame import convert
from narupatools.lammps import LAMMPSSimulation

from ._calculator import LAMMPSCalculator


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


__all__ = ["atoms_from_lammps_simulation", "LAMMPSCalculator"]
