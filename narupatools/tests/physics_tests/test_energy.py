import pytest
from numpy.random import default_rng

from narupatools.physics.energy import cumulative_work, total_work, work_per_step

rng = default_rng()


@pytest.mark.parametrize(
    ("func",), [(total_work,), (cumulative_work,), (work_per_step,)]
)
def test_work_incompatible_components(func):
    forces = rng.standard_normal((10, 2))
    positions = rng.standard_normal((10, 3))
    with pytest.raises(ValueError):
        _ = func(forces=forces, positions=positions)


@pytest.mark.parametrize(
    ("func",), [(total_work,), (cumulative_work,), (work_per_step,)]
)
def test_work_incompatible_timestep(func):
    forces = rng.standard_normal((10, 3))
    positions = rng.standard_normal((15, 3))
    with pytest.raises(ValueError):
        _ = func(forces=forces, positions=positions)


@pytest.mark.parametrize(
    ("func",), [(total_work,), (cumulative_work,), (work_per_step,)]
)
def test_work_incompatible_grouping(func):
    forces = rng.standard_normal((15, 10, 2))
    positions = rng.standard_normal((5, 10, 3))
    with pytest.raises(ValueError):
        _ = func(forces=forces, positions=positions)


@pytest.mark.parametrize(
    ("func", "shape"),
    [(total_work, ()), (cumulative_work, (10,)), (work_per_step, (10,))],
)
def test_work_shape(func, shape):
    forces = rng.standard_normal((10, 3))
    positions = rng.standard_normal((10, 3))
    work = func(forces=forces, positions=positions)
    assert work.shape == shape


@pytest.mark.parametrize(
    ("func", "shape"),
    [(total_work, (5,)), (cumulative_work, (10, 5)), (work_per_step, (10, 5))],
)
def test_work_multiple_particles_shape(func, shape):
    forces = rng.standard_normal((10, 5, 3))
    positions = rng.standard_normal((10, 5, 3))
    work = func(forces=forces, positions=positions, time_axis=-3)
    assert work.shape == shape


@pytest.mark.parametrize(
    ("func", "shape"),
    [(total_work, (15,)), (cumulative_work, (15, 10)), (work_per_step, (15, 10))],
)
def test_total_work_grouping_shape(func, shape):
    forces = rng.standard_normal((15, 10, 3))
    positions = rng.standard_normal((10, 3))
    work = func(forces=forces, positions=positions)
    assert work.shape == shape


@pytest.mark.parametrize(
    ("func", "shape"),
    [
        (total_work, (15, 5)),
        (cumulative_work, (15, 10, 5)),
        (work_per_step, (15, 10, 5)),
    ],
)
def test_total_work_multiple_particle_grouping_shape(func, shape):
    forces = rng.standard_normal((15, 10, 5, 3))
    positions = rng.standard_normal((10, 5, 3))
    work = func(forces=forces, positions=positions, time_axis=-3)
    assert work.shape == shape
