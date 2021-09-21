import random

import numpy as np
import pytest
from _pytest.python_api import ApproxBase

from narupatools.ase import ASEDynamics, ConstantCalculator, NullCalculator, UnitsASE
from narupatools.ase._rotational_velocity_verlet import (
    RotationalVelocityVerletIntegrator,
)
from narupatools.ase._system import ASESystem
from narupatools.core import UnitsNarupa
from narupatools.core.random import random_float, random_integer
from narupatools.imd import rigidmotion_interaction
from narupatools.physics._quaternion import from_rotation_vector, quaternion
from narupatools.physics.random import (
    random_unit_quaternion,
    random_vector,
)
from narupatools.physics.transformation import Rotation
from narupatools.physics.vector import normalized, vector

_NarupaToASE = UnitsNarupa >> UnitsASE


def test_stationary(single_carbon_atoms):
    atoms = single_carbon_atoms
    atoms.calc = NullCalculator()
    integrator = RotationalVelocityVerletIntegrator(
        atoms, timestep=0.1 * _NarupaToASE.time
    )

    dynamics = ASEDynamics(integrator)

    initial = dynamics.positions

    integrator.run(100)

    assert dynamics.positions == pytest.approx(initial)
    assert dynamics.velocities == pytest.approx(np.array([[0.0, 0.0, 0.0]]))


@pytest.fixture(params=range(1))
def seed(request):
    random.seed(request.param)
    return request.param


@pytest.fixture
def velocity(seed):
    return random_vector()


@pytest.fixture
def timestep(seed):
    return random_float(minimum=0.001, maximum=0.01)


@pytest.fixture
def nsteps(seed):
    return random_integer(minimum=10, maximum=500)


@pytest.fixture
def mass(seed, single_carbon_atoms):
    mass = np.array([random_float(minimum=0.01, maximum=1.0)])
    ASESystem(single_carbon_atoms).masses = mass
    return mass


@pytest.fixture
def symmetric_inertia(seed, single_carbon_atoms):
    inertia = np.array([random_float(minimum=0.01, maximum=1.0)])
    ASESystem(single_carbon_atoms).moments_of_inertia = inertia
    return inertia


@pytest.fixture
def angular_velocities(seed, single_carbon_atoms):
    omega = np.array([random_vector()])
    ASESystem(single_carbon_atoms).angular_velocities = omega
    return omega


@pytest.fixture
def orientations(seed, single_carbon_atoms):
    orientations = np.array([random_unit_quaternion()], dtype=quaternion)
    ASESystem(single_carbon_atoms).orientations = orientations
    return orientations


class ApproxQuaternionArray(ApproxBase):
    def __repr__(self) -> str:
        return f"approx({self.expected})"

    def __eq__(self, actual) -> bool:
        if actual.dtype != quaternion:
            return False

        if actual.shape != self.expected.shape:
            return False

        return ApproxBase.__eq__(self, actual)

    def _yield_comparisons(self, actual):
        for i in range(self.expected.shape[0]):
            quat_actual = actual[i]
            quat_expected = self.expected[i]
            for i in range(4):
                yield quat_actual.components[i], quat_expected.components[i]


def approx(obj, rel=None, abs=None):  # noqa: A002
    if obj.dtype == quaternion:
        return ApproxQuaternionArray(obj, rel, abs)
    else:
        return pytest.approx(obj)


def test_no_forces(
        mass,
        symmetric_inertia,
        angular_velocities,
        orientations,
        velocity,
        timestep,
        nsteps,
        single_carbon_atoms,
):
    atoms = single_carbon_atoms
    atoms.calc = NullCalculator()
    integrator = RotationalVelocityVerletIntegrator(
        atoms, timestep=timestep * _NarupaToASE.time
    )

    dynamics = ASEDynamics(integrator)
    dynamics.velocities = np.array([velocity])

    initial = dynamics.positions

    integrator.run(nsteps)

    assert dynamics.positions == approx(initial + velocity * nsteps * timestep)

    assert dynamics.orientations == approx(
        from_rotation_vector(angular_velocities * nsteps * timestep) * orientations,
        rel=1e-3,
        abs=1e-3,
    )


def test_constant_torque(
        mass, symmetric_inertia, angular_velocities, nsteps, timestep, single_carbon_atoms
):
    atoms = single_carbon_atoms
    torque = random_vector(max_magnitude=0.5)

    atoms.calc = ConstantCalculator(
        forces=[vector(0, 0, 0)], torques=[torque * _NarupaToASE.torque]
    )

    integrator = RotationalVelocityVerletIntegrator(
        atoms, timestep=timestep * _NarupaToASE.time
    )

    dynamics = ASEDynamics(integrator)

    start = dynamics.angular_momenta

    start_angvel = dynamics.angular_velocities

    dynamics.run(nsteps)

    assert dynamics.angular_momenta == pytest.approx(start + nsteps * timestep * torque)
    assert dynamics.angular_velocities == pytest.approx(
        start_angvel + nsteps * timestep * torque / dynamics.moments_of_inertia[0]
    )


def test_rotate(mass, symmetric_inertia, timestep, single_carbon_atoms):
    atoms = single_carbon_atoms
    atoms.calc = NullCalculator()
    integrator = RotationalVelocityVerletIntegrator(
        atoms, timestep=timestep * _NarupaToASE.time
    )

    dynamics = ASEDynamics(integrator)

    dynamics.moments_of_inertia = [2.0]

    rotation = normalized(quaternion(1, 1, 0, 0))

    dynamics.imd.add_interaction(
        rigidmotion_interaction(particles=[0], scale=100.0, rotation=rotation)
    )

    for _ in range(100):
        print(Rotation(dynamics.orientations[0]))
        print(dynamics.torques[0])
        dynamics.run(5)

    end = dynamics.orientations[0]

    print(rotation)
    print(end)
