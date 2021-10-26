import random

import numpy as np
import pytest
from ase import Atoms

from narupatools.ase import ASEDynamics
from narupatools.imd import rigidmotion_interaction
from narupatools.physics.transformation import Translation
from narupatools.physics.vector import vector


@pytest.fixture
def methane_positions():
    return (
        np.array(
            [
                [0, 0, 0],
                [0, 0, 1.089],
                [1.027, 0, -0.363],
                [-0.513, -0.889, -0.363],
                [-0.513, 0.889, -0.363],
            ]
        )
        * 0.1
    )


@pytest.fixture
def single_carbon_atoms(methane_positions) -> Atoms:
    return Atoms(
        symbols=["C", "H", "H", "H", "H"],
        masses=np.array([12.0, 1.0, 1.0, 1.0, 1.0]),
        positions=methane_positions * 10.0,
    )


@pytest.fixture
def dynamics(single_carbon_atoms) -> ASEDynamics:
    return ASEDynamics.create_velocity_verlet(single_carbon_atoms, timestep=0.1)


@pytest.fixture(params=range(20))
def seed(request):
    random.seed(request.param)
    return request.param


def test_translate_zero(dynamics, methane_positions):
    interaction = rigidmotion_interaction(
        particles=[0, 1, 2, 3, 4], translation=vector(0, 0, 0)
    )
    dynamics.imd.add_interaction(interaction)
    dynamics.run(100)
    assert dynamics.positions == pytest.approx(methane_positions)


@pytest.fixture(
    params=[
        vector(0, 0, 0),
        vector(2, 0, 0),
        vector(0, 5, 0),
        vector(0, 0, -8),
        vector(2, -3, -1),
    ]
)
def translation(request):
    return request.param


def test_translate_1(dynamics, methane_positions, translation):
    interaction = rigidmotion_interaction(
        particles=[0, 1, 2, 3, 4], translation=translation, scale=100
    )
    dynamics.imd.add_interaction(interaction)
    dynamics.run(500)
    assert dynamics.positions == pytest.approx(
        Translation(translation) @ methane_positions, rel=1e-1
    )
