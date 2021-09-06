import random

import numpy as np
import pytest
from _pytest.python_api import _recursive_list_map, ApproxBase

from narupatools.ase import UnitsASE, NullCalculator, ASEDynamics
from narupatools.ase._rotational_velocity_verlet import RotationalVelocityVerletIntegrator
from narupatools.ase._system import ASESystem
from narupatools.core import UnitsNarupa
from narupatools.core.random import random_integer, random_float
from narupatools.imd import rigidmotion_interaction
from narupatools.physics.random import random_vector, random_quaternion, random_unit_quaternion
from narupatools.physics.quaternion import from_rotation_vector, quaternion
from narupatools.physics.transformation import Rotation

_NarupaToASE = UnitsNarupa >> UnitsASE


def test_stationary(single_carbon_atoms):
    atoms = single_carbon_atoms
    atoms.calc = NullCalculator()
    integrator = RotationalVelocityVerletIntegrator(atoms, timestep=0.1 * _NarupaToASE.time)

    dynamics = ASEDynamics(integrator)

    initial = dynamics.positions

    integrator.run(100)

    assert dynamics.positions == pytest.approx(initial)
    assert dynamics.velocities == pytest.approx(np.array([[0.0, 0.0, 0.0]]))


@pytest.fixture(params=range(20))
def seed(request):
    random.seed(request.param)
    return request.param


@pytest.fixture
def velocity(seed):
    return random_vector()


@pytest.fixture
def timestep(seed):
    return random_float(min=0.001, max=0.01)


@pytest.fixture
def nsteps(seed):
    return random_integer(min=10, max=500)


@pytest.fixture
def mass(seed, single_carbon_atoms):
    mass = np.array([random_float(min=0.01, max=1.0)])
    ASESystem(single_carbon_atoms).masses = mass
    return mass


@pytest.fixture
def symmetric_inertia(seed, single_carbon_atoms):
    inertia = np.array([random_float(min=0.01, max=1.0)])
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


def approx(obj, rel = None, abs = None):
    if obj.dtype == quaternion:
        return ApproxQuaternionArray(obj, rel, abs)
    else:
        return pytest.approx(obj)


def test_no_forces(mass, symmetric_inertia, angular_velocities, orientations, velocity, timestep, nsteps, single_carbon_atoms):
    atoms = single_carbon_atoms
    atoms.calc = NullCalculator()
    integrator = RotationalVelocityVerletIntegrator(atoms, timestep=timestep * _NarupaToASE.time)

    dynamics = ASEDynamics(integrator)
    dynamics.velocities = np.array([velocity])

    initial = dynamics.positions

    integrator.run(nsteps)

    assert dynamics.positions == approx(initial + velocity * nsteps * timestep)

    assert dynamics.orientations == approx(from_rotation_vector(angular_velocities * nsteps * timestep) * orientations, rel=1e-3, abs=1e-3)



def test_rotate(mass, symmetric_inertia, timestep, single_carbon_atoms):
    atoms = single_carbon_atoms
    atoms.calc = NullCalculator()
    integrator = RotationalVelocityVerletIntegrator(atoms, timestep = timestep * _NarupaToASE.time)

    dynamics = ASEDynamics(integrator)

    rotation = random_unit_quaternion()

    start = dynamics.orientations[0]

    dynamics.imd.add_interaction(rigidmotion_interaction(particles=[0], scale=10.0, rotation=rotation))

    dynamics.run(100)

    end = dynamics.orientations[0]