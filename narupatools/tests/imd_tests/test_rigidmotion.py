from narupatools.imd.interactions import rigidmotion_interaction


def test_no_translation_rotation():
    data = rigidmotion_interaction(particles=[0, 3])

    assert not hasattr(data, "rotation")
    assert not hasattr(data, "position")
