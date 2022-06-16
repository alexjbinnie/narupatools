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

import math

import numpy as np
import pytest

from narupatools.physics.force import (
    centripetal_force,
    critically_damped_spring_force,
    damped_spring_force,
    gaussian_force_and_energy,
    linear_drag_force,
    mass_weighted_forces,
    spring_force,
    spring_force_and_energy,
)
from narupatools.physics.random import random_float, random_integer, random_vector
from narupatools.physics.vector import (
    dot_product,
    magnitude,
    sqr_magnitude,
    vector_rejection,
    zero_vector,
)


@pytest.fixture
def offset(seed):
    return random_vector(max_magnitude=100.0)


@pytest.fixture
def spring_constant(seed):
    return random_float(minimum=0.0, maximum=100.0)


@pytest.fixture
def depth(seed):
    return random_float(minimum=0.0, maximum=100.0)


@pytest.fixture
def sigma(seed):
    return random_float(minimum=0.0, maximum=100.0)


@pytest.fixture
def velocity(seed):
    return random_vector(max_magnitude=100.0)


@pytest.fixture
def position(seed):
    return random_vector(max_magnitude=100.0)


@pytest.fixture
def origin(seed):
    return random_vector(max_magnitude=100.0)


@pytest.fixture
def angular_velocity(seed):
    return random_vector(max_magnitude=100.0)


@pytest.fixture
def force(seed):
    return random_vector(max_magnitude=100.0)


@pytest.fixture
def damping_coefficient(seed):
    return random_float(minimum=0.0, maximum=100.0)


@pytest.fixture
def mass(seed):
    return random_float(minimum=0.0, maximum=100.0)


@pytest.fixture
def system_size(seed):
    return random_integer(minimum=1, maximum=10)


@pytest.fixture
def masses(seed, system_size):
    return [random_float(minimum=0.0, maximum=100.0) for _ in range(system_size)]


@pytest.fixture
def positions(seed, system_size):
    return [random_vector(max_magnitude=100.0) for _ in range(system_size)]


def test_spring_force_zero_offset(spring_constant):
    assert spring_force(
        offset=zero_vector(), spring_constant=spring_constant
    ) == pytest.approx(zero_vector())


def test_spring_force_magnitude(offset, spring_constant):
    assert magnitude(
        spring_force(offset=offset, spring_constant=spring_constant)
    ) == pytest.approx(spring_constant * magnitude(offset))


def test_spring_force_opposes_offset(offset):
    assert dot_product(spring_force(offset=offset, spring_constant=1.0), offset) < 0


def test_spring_force_and_energy_same_force(offset, spring_constant):
    force, _ = spring_force_and_energy(offset=offset, spring_constant=spring_constant)
    assert spring_force(
        offset=offset, spring_constant=spring_constant
    ) == pytest.approx(force)


def test_spring_energy_zero_offset(spring_constant):
    _, energy = spring_force_and_energy(
        offset=zero_vector(), spring_constant=spring_constant
    )
    assert energy == pytest.approx(0.0)


def test_spring_energy_correct(offset, spring_constant):
    _, energy = spring_force_and_energy(offset=offset, spring_constant=spring_constant)
    assert energy == pytest.approx(0.5 * (magnitude(offset) ** 2.0) * spring_constant)


def test_gaussian_energy_zero_offset(depth, sigma):
    _, energy = gaussian_force_and_energy(
        offset=zero_vector(), depth=depth, sigma=sigma
    )
    assert energy == pytest.approx(-depth)


def test_gaussian_force_zero_offset(depth, sigma):
    force, _ = gaussian_force_and_energy(offset=zero_vector(), depth=depth, sigma=sigma)
    assert force == pytest.approx(zero_vector())


def test_gaussian_energy(offset, depth, sigma):
    _, energy = gaussian_force_and_energy(offset=offset, depth=depth, sigma=sigma)
    assert energy == pytest.approx(
        -depth * math.exp(-sqr_magnitude(offset) / (2 * sigma**2))
    )


def test_gaussian_force(offset, depth, sigma):
    force, _ = gaussian_force_and_energy(offset=offset, depth=depth, sigma=sigma)
    assert force == pytest.approx(
        -depth
        * offset
        * math.exp(-sqr_magnitude(offset) / (2 * sigma**2))
        / (sigma**2)
    )


def test_linear_drag_force(velocity, damping_coefficient):
    force = linear_drag_force(
        velocity=velocity, damping_coefficient=damping_coefficient
    )
    assert force == pytest.approx(-damping_coefficient * velocity)


def test_damped_spring_force(*, offset, velocity, spring_constant, damping_coefficient):
    force = damped_spring_force(
        offset=offset,
        velocity=velocity,
        spring_constant=spring_constant,
        damping_coefficient=damping_coefficient,
    )
    assert force == pytest.approx(
        -spring_constant * offset - damping_coefficient * velocity
    )


def test_critically_damped_spring_force(*, offset, velocity, spring_constant, mass):
    force = critically_damped_spring_force(
        offset=offset, velocity=velocity, spring_constant=spring_constant, mass=mass
    )
    damping_coefficient = 2 * math.sqrt(mass * spring_constant)
    assert force == pytest.approx(
        -spring_constant * offset - damping_coefficient * velocity
    )


def test_mass_weighted_forces(masses, force):
    forces = mass_weighted_forces(masses=masses, force=force)
    total_mass = sum(masses)
    for i, mass in enumerate(masses):
        assert forces[i] == pytest.approx(mass * force / total_mass)


def test_mass_weighted_forces_equal_accelerations(masses, force):
    forces = mass_weighted_forces(masses=masses, force=force)
    total_mass = sum(masses)
    acceleration = force / total_mass
    for i, mass in enumerate(masses):
        if mass > 0:
            assert forces[i] / mass == pytest.approx(acceleration)


def test_centripetal_force_single_particle_is_zero(mass, position, angular_velocity):
    force = centripetal_force(
        masses=[mass], positions=[position], angular_velocity=angular_velocity
    )
    assert force == pytest.approx(np.array([[0.0, 0.0, 0.0]]), abs=1e-6)

    force2 = centripetal_force(
        masses=[mass],
        positions=[position],
        angular_velocity=angular_velocity,
        origin=position,
    )
    assert force2 == pytest.approx(np.array([[0.0, 0.0, 0.0]]), abs=1e-6)


def test_centripetal_force_single_particle_arbitrary_axis(
    mass, position, angular_velocity, origin
):
    force = centripetal_force(
        masses=[mass],
        positions=[position],
        angular_velocity=angular_velocity,
        origin=origin,
    )
    rperp = vector_rejection(np.subtract(position, origin), angular_velocity)
    assert force == pytest.approx(
        np.asfarray([-mass * sqr_magnitude(angular_velocity) * rperp])
    )
