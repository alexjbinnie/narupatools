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
