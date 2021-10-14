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

import pytest

lammps = pytest.importorskip("lammps")

from narupatools.ase import ASEDynamics
from narupatools.ase.lammps import atoms_from_lammps_simulation
from narupatools.lammps import LAMMPSDynamics, LAMMPSSimulation


def test_ase_lammps():
    simulation = LAMMPSSimulation.from_file("in.peptide")

    atoms = atoms_from_lammps_simulation(simulation)

    dynamics = ASEDynamics.create_velocity_verlet(atoms, timestep=0.002)

    dynamics.run(100)


def test_lammps():
    simulation = LAMMPSSimulation.from_file("in.peptide")

    dynamics = LAMMPSDynamics(simulation)

    dynamics.run(100)
