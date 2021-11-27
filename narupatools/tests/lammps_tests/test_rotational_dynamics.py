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

import random

import numpy as np
import pytest

lammps = pytest.importorskip("lammps")

from narupatools.lammps import INSTALLED_PACKAGES, LAMMPSSimulation
from narupatools.lammps.regions import Box
from narupatools.physics.random import random_float, random_integer, random_vector
from narupatools.physics.transformation import Rotation
from narupatools.physics.vector import vector

if "ASPHERE" not in INSTALLED_PACKAGES:
    pytest.skip("Requires ASPHERE package.", allow_module_level=True)

if "CG-DNA" not in INSTALLED_PACKAGES:
    pytest.skip("Requires CG-DNA package.", allow_module_level=True)

@pytest.fixture(params=range(20))
def seed(request):
    random.seed(request.param)
    return request.param


@pytest.fixture
def radius(seed):
    return random_float(minimum=0.2, maximum=10.0)


@pytest.fixture
def timestep(seed):
    return random_float(minimum=0.01, maximum=0.1)


@pytest.fixture
def simulation(timestep):
    simulation = LAMMPSSimulation.create_new("nano")
    simulation.command("atom_style ellipsoid")
    simulation.create_box(n_types=1, region=Box.from_size(size=vector(10, 10, 10)))
    simulation.timestep = timestep
    simulation.command("fix integrate all nve/dot")
    return simulation


@pytest.fixture
def angular_momentum(seed):
    return random_vector(max_magnitude=10.0)


@pytest.fixture
def mass(seed):
    return random_float(minimum=0.1, maximum=10.0)


@pytest.fixture
def spherical_atom(seed, simulation, radius, mass, angular_momentum):
    atom = simulation.create_atom(
        atom_type=1, position=vector(0, 0, 0), rotation=Rotation.identity
    )
    atom.set_shape(shape=vector(1, 1, 1) * radius * 2.0)
    atom.set_peratom_mass(mass=mass)
    return atom


@pytest.fixture
def initial_angular_momentum(spherical_atom, angular_momentum):
    spherical_atom.set_angular_momentum(angular_momentum=angular_momentum)


@pytest.fixture
def n_steps(seed):
    return random_integer(minimum=5, maximum=50)


def test_spherical_angmom(
    simulation,
    mass,
    radius,
    angular_momentum,
    spherical_atom,
    initial_angular_momentum,
    n_steps,
):
    simulation.run(n_steps)

    # Moment of inertia of a sphere
    moment_inertia = 0.4 * mass * (radius ** 2)

    angular_displacement = (
        n_steps * simulation.timestep * angular_momentum / moment_inertia
    )

    assert Rotation.from_rotation_vector(
        angular_displacement
    ).versor.components == pytest.approx(
        Rotation(simulation.orientations[0]).versor.components, rel=1e-2, abs=1e-1
    )


def test_spherical_torque(simulation, mass, radius, spherical_atom, n_steps):
    torque = random_vector(max_magnitude=10.0)

    def callback():
        return 0, None, np.array([torque])

    simulation.setup_imd(callback)
    simulation.run(0, pre=True)

    simulation.run(n_steps)

    angmom_calc = (n_steps * simulation.timestep) * torque
    angmom_act = simulation.angular_momenta[0]

    assert angmom_act == pytest.approx(angmom_calc)
