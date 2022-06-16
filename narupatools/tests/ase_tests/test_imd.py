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

import numpy as np
import pytest
from ase.atoms import Atoms
from ase.md import VelocityVerlet

from narupatools.ase import ASEDynamics, NullCalculator, UnitsASE
from narupatools.imd import gaussian_interaction
from narupatools.physics.units import UnitsNarupa

_NarupaToASE = UnitsNarupa >> UnitsASE


@pytest.fixture
def single_carbon_atoms():
    atoms = Atoms(symbols=["C"], positions=[[0, 0, 0]])
    atoms.set_calculator(NullCalculator())
    return atoms


@pytest.fixture
def single_carbon_verlet(single_carbon_atoms):
    return ASEDynamics(VelocityVerlet(single_carbon_atoms, timestep=0.01))


@pytest.fixture
def single_carbon_imd(single_carbon_verlet):
    imd = {}
    single_carbon_verlet.imd.add_source(imd)
    return imd


def test_system_stationary(single_carbon_atoms, single_carbon_verlet):
    single_carbon_verlet.run(steps=50)
    assert np.linalg.norm(single_carbon_atoms.get_velocities()[0]) == pytest.approx(0.0)


def test_initial_no_interactions(single_carbon_verlet):
    assert len(single_carbon_verlet.imd.current_interactions) == 0


def test_add_interaction(single_carbon_imd, single_carbon_verlet):
    single_carbon_imd["abc"] = gaussian_interaction(particles=[0], position=[1, 0, 0])
    assert len(single_carbon_verlet.imd.current_interactions) == 0
    single_carbon_verlet.run(steps=1)
    assert len(single_carbon_verlet.imd.current_interactions) == 1


def test_interaction_applies_force(
    single_carbon_imd,
    single_carbon_atoms,
    single_carbon_verlet,
):
    single_carbon_imd["abc"] = gaussian_interaction(
        particles=[0], position=[1, 0, 0], scale=1000.0
    )
    single_carbon_verlet.run(steps=50)
    assert abs(single_carbon_atoms.get_velocities()[0][0]) > 0.0
    assert abs(single_carbon_atoms.get_velocities()[0][1]) == pytest.approx(0.0)
    assert abs(single_carbon_atoms.get_velocities()[0][2]) == pytest.approx(0.0)


def test_work_is_kinetic_energy(
    single_carbon_imd,
    single_carbon_atoms,
    single_carbon_verlet,
):
    single_carbon_imd["abc"] = gaussian_interaction(
        particles=[0], position=[1, 0, 0], scale=100.0
    )
    assert single_carbon_atoms.get_kinetic_energy() == pytest.approx(0.0)
    single_carbon_verlet.run(steps=100)
    work_done = single_carbon_verlet.imd.total_work
    assert single_carbon_atoms.get_kinetic_energy() == pytest.approx(
        work_done * _NarupaToASE.energy, rel=1e-2
    )


def test_interaction_pickle(
    single_carbon_imd,
    single_carbon_atoms,
    single_carbon_verlet,
):
    single_carbon_imd["abc"] = gaussian_interaction(
        particles=[0], position=[1, 0, 0], scale=1000.0
    )
    single_carbon_verlet.run(1)
    _ = single_carbon_atoms.copy()
