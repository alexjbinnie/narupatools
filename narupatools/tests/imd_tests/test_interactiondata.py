import numpy as np
import pytest
from narupa.imd import ParticleInteraction
from narupa.imd.imd_state import dict_to_interaction, interaction_to_dict

from narupatools.imd.interactions import GAUSSIAN_INTERACTION_TYPE
from narupatools.imd.interactions._interactiondata import InteractionData
from narupatools.imd.interactions._point import (
    SPRING_INTERACTION_TYPE,
    PointInteractionData,
)


def test_gaussian_particle_interaction_to_interactiondata():
    particle_interaction = ParticleInteraction(
        interaction_type=GAUSSIAN_INTERACTION_TYPE,
        particles=[0, 3, 5],
        scale=2.5,
        position=[0.2, 1.0, -2.0],
    )

    dict = interaction_to_dict(particle_interaction)

    data = InteractionData.deserialize(dict)

    assert data.interaction_type == GAUSSIAN_INTERACTION_TYPE
    assert isinstance(data, PointInteractionData)
    assert data.particles == pytest.approx(np.array([0, 3, 5]))
    assert data.scale == pytest.approx(2.5)
    assert data.position == pytest.approx([0.2, 1.0, -2.0])

    dict2 = data.serialize()

    particle_interaction2 = dict_to_interaction(dict2)

    assert particle_interaction == particle_interaction2


def test_spring_particle_interaction_to_interactiondata():
    particle_interaction = ParticleInteraction(
        interaction_type=SPRING_INTERACTION_TYPE,
        particles=[0, 3, 5],
        scale=2.5,
        position=[0.2, 1.0, -2.0],
    )

    dict = interaction_to_dict(particle_interaction)

    data = InteractionData.deserialize(dict)

    assert data.interaction_type == SPRING_INTERACTION_TYPE
    assert isinstance(data, PointInteractionData)
    assert data.particles == pytest.approx(np.array([0, 3, 5]))
    assert data.scale == pytest.approx(2.5)
    assert data.position == pytest.approx([0.2, 1.0, -2.0])

    dict2 = data.serialize()

    particle_interaction2 = dict_to_interaction(dict2)

    assert particle_interaction == particle_interaction2
