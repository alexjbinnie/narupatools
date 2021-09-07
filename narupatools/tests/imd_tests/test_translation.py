import math

import numpy as np
import pytest
from ase import Atoms

from narupatools.ase import ASEDynamics
from narupatools.frame.hdf5 import add_hdf5_writer
from narupatools.imd import rigidmotion_interaction
from narupatools.physics.transformation import Rotation, Translation
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
        Translation(translation) * methane_positions
    )


def test_rotation(dynamics, methane_positions):
    angle = vector(0.5 * math.pi, 0, 0)

    interaction = rigidmotion_interaction(
        particles=[0, 1, 2, 3, 4], rotation=angle, scale=15
    )
    dynamics.imd.add_interaction(interaction)

    dynamics.run(500)

    rotation = Rotation.from_rotation_vector(angle)

    np.set_printoptions(suppress=True)

    assert dynamics.positions == pytest.approx(
        rotation * methane_positions, rel=1e-2, abs=1e-2
    )


def test_rotation_translation(dynamics, methane_positions):
    writer = add_hdf5_writer(dynamics, filename="test.hdf5", overwrite_existing=True)

    angle = vector(0.5 * math.pi, 0, 0)
    translate = vector(0.4, 0, 0)

    interaction = rigidmotion_interaction(
        particles=[0, 1, 2, 3, 4], translation=translate, rotation=angle, scale=5
    )

    dynamics.imd.add_interaction(interaction)

    rotation = Rotation.from_rotation_vector(angle)
    translation = Translation(translate)

    dynamics.run(200)

    assert dynamics.positions == pytest.approx(
        translation * (rotation * methane_positions), rel=1e-2, abs=1e-2
    )

    writer.close()
