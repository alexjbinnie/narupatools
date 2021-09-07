import numpy as np
import pytest

from narupatools.core.random import random_integer
from narupatools.physics import rigidbody
from narupatools.physics.matrix import identity_matrix, kronecker_delta, zero_matrix
from narupatools.physics.random import random_float, random_quaternion, random_vector
from narupatools.physics.rigidbody import (
    center_of_mass,
    center_of_mass_acceleration,
    center_of_mass_velocity,
    distribute_angular_velocity,
    kinetic_energy,
    moment_of_inertia_tensor,
    orbital_angular_momentum,
    spin_angular_momentum,
)
from narupatools.physics.vector import (
    cross_product,
    dot_product,
    outer_product,
    sqr_magnitude,
    zero_vector,
)


@pytest.fixture
def system_size(seed):
    return random_integer(1, 10)


@pytest.fixture
def positions(seed, system_size):
    return [random_vector(max_magnitude=100.0) for _ in range(system_size)]


@pytest.fixture
def velocities(seed, system_size):
    return [random_vector(max_magnitude=100.0) for _ in range(system_size)]


@pytest.fixture
def accelerations(seed, system_size):
    return [random_vector(max_magnitude=100.0) for _ in range(system_size)]


@pytest.fixture
def masses(seed, system_size):
    return [random_float(max=100.0) for _ in range(system_size)]


@pytest.fixture
def origin(seed):
    return random_vector(max_magnitude=100.0)


@pytest.fixture
def displacement(seed):
    return random_vector(max_magnitude=100.0)


@pytest.fixture
def rotation(seed):
    return [random_quaternion(max_magnitude=100.0) for _ in range(system_size)]


@pytest.fixture
def angular_velocity(seed):
    return random_vector(max_magnitude=100.0)


def test_center_of_mass(masses, positions):
    com = center_of_mass(positions=positions, masses=masses)

    assert com.shape == (3,)

    com_calc = zero_vector()
    for position, mass in zip(positions, masses):
        com_calc += np.asfarray(position) * mass
    com_calc /= sum(masses)

    assert com == pytest.approx(com_calc)


def test_center_of_mass_velocity(masses, velocities):
    com = center_of_mass_velocity(velocities=velocities, masses=masses)

    assert com.shape == (3,)

    com_calc = zero_vector()
    for velocity, mass in zip(velocities, masses):
        com_calc += np.asfarray(velocity) * mass
    com_calc /= sum(masses)

    assert com == pytest.approx(com_calc)


def test_center_of_mass_acceleration(masses, accelerations):
    com = center_of_mass_acceleration(accelerations=accelerations, masses=masses)

    assert com.shape == (3,)

    com_calc = zero_vector()
    for acceleration, mass in zip(accelerations, masses):
        com_calc += np.asfarray(acceleration) * mass
    com_calc /= sum(masses)

    assert com == pytest.approx(com_calc)


def test_spin_angular_momentum(masses, positions, velocities):
    angular_momentum = spin_angular_momentum(
        masses=masses, positions=positions, velocities=velocities
    )

    assert angular_momentum.shape == (3,)

    com = center_of_mass(masses=masses, positions=positions)
    com_vel = center_of_mass_velocity(masses=masses, velocities=velocities)

    ang_mom_calc = zero_vector()
    for mass, position, velocity in zip(masses, positions, velocities):
        ang_mom_calc += mass * cross_product(position - com, velocity - com_vel)

    assert angular_momentum == pytest.approx(ang_mom_calc)


def test_orbit_angular_momentum(masses, positions, velocities):
    angular_momentum = orbital_angular_momentum(
        masses=masses, positions=positions, velocities=velocities
    )

    assert angular_momentum.shape == (3,)

    com = center_of_mass(masses=masses, positions=positions)
    com_vel = center_of_mass_velocity(masses=masses, velocities=velocities)

    ang_mom_calc = sum(masses) * cross_product(com, com_vel)

    assert angular_momentum == pytest.approx(ang_mom_calc)


def test_orbit_angular_momentum_origin(masses, positions, velocities, origin):
    angular_momentum = orbital_angular_momentum(
        masses=masses, positions=positions, velocities=velocities, origin=origin
    )

    assert angular_momentum.shape == (3,)

    com = center_of_mass(masses=masses, positions=positions)
    com_vel = center_of_mass_velocity(masses=masses, velocities=velocities)

    ang_mom_calc = sum(masses) * cross_product(com - origin, com_vel)

    assert angular_momentum == pytest.approx(ang_mom_calc)


def test_total_angular_momentum(masses, positions, velocities):
    total = spin_angular_momentum(
        masses=masses, positions=positions, velocities=velocities
    ) + orbital_angular_momentum(
        masses=masses, positions=positions, velocities=velocities
    )

    calc = zero_vector()
    for mass, position, velocity in zip(masses, positions, velocities):
        calc += mass * cross_product(position, velocity)

    assert total == pytest.approx(calc)


def test_inertia_tensor(masses, positions):
    inertia_tensor = moment_of_inertia_tensor(masses=masses, positions=positions)
    tensor_calc = zero_matrix()

    com = center_of_mass(masses=masses, positions=positions)
    for i in range(3):
        for j in range(3):
            for mass, position in zip(masses, positions):
                offset = position - com
                tensor_calc[i][j] += mass * (
                    sqr_magnitude(offset) * kronecker_delta(i, j)
                    - offset[i] * offset[j]
                )

    assert inertia_tensor == pytest.approx(tensor_calc)


def test_inertia_tensor_translation(masses, positions, displacement):
    inertia_tensor = moment_of_inertia_tensor(masses=masses, positions=positions)
    com = center_of_mass(masses=masses, positions=positions)
    inertia_tensor_displaced = moment_of_inertia_tensor(
        masses=masses, positions=positions, origin=com + displacement
    )
    total_mass = sum(masses)
    assert inertia_tensor_displaced == pytest.approx(
        inertia_tensor
        + total_mass
        * (
            dot_product(displacement, displacement) * identity_matrix()
            - outer_product(displacement, displacement)
        )
    )


def test_distribute_angular_velocity_and_reverse(
    system_size, masses, positions, angular_velocity
):
    if system_size <= 2:
        return
    velocities = distribute_angular_velocity(
        angular_velocity=angular_velocity, masses=masses, positions=positions
    )
    calc_angular_velocity = rigidbody.angular_velocity(
        masses=masses, positions=positions, velocities=velocities
    )
    assert calc_angular_velocity == pytest.approx(angular_velocity)


def test_distribute_angular_velocity(masses, positions, angular_velocity):
    velocities = distribute_angular_velocity(
        angular_velocity=angular_velocity, masses=masses, positions=positions
    )
    com = center_of_mass(masses=masses, positions=positions)

    for i, position in enumerate(positions):
        assert velocities[i] == pytest.approx(
            cross_product(angular_velocity, position - com)
        )


def test_distribute_angular_velocity_origin(
    masses, positions, angular_velocity, origin
):
    velocities = distribute_angular_velocity(
        angular_velocity=angular_velocity,
        masses=masses,
        positions=positions,
        origin=origin,
    )

    for i, position in enumerate(positions):
        assert velocities[i] == pytest.approx(
            cross_product(angular_velocity, position - origin)
        )


def test_distribute_angular_velocity_no_origin_or_mass(positions, angular_velocity):
    with pytest.raises(ValueError):  # noqa: PT011
        _ = distribute_angular_velocity(
            angular_velocity=angular_velocity, positions=positions
        )


def test_kinetic_energy(masses, velocities):
    energy = kinetic_energy(masses=masses, velocities=velocities)

    energy_calc = 0.0
    for mass, velocity in zip(masses, velocities):
        energy_calc += 0.5 * mass * sqr_magnitude(velocity)

    assert energy == pytest.approx(energy_calc)
