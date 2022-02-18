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

from narupatools.physics.integrate import (
    cumulative_vector_line_integral,
    integral,
    vector_line_integral,
)


def samples(start, end, sample_count):
    return start + (np.arange(sample_count) / (sample_count - 1)) * (end - start)


def test_samples():
    assert samples(2, 4, 5) == pytest.approx([2, 2.5, 3, 3.5, 4])


INTEGRALS = [
    (lambda x: x**2, 0, 1, 100, 0.33335034),
    (lambda x: x**2, 2, 4, 5, 18.75),
    (lambda x: math.sin(x), 0, math.pi, 50, 1.9993148),
    (lambda x: math.exp(x), -3, 3, 80, 20.04538),
    (lambda x: math.cosh(x), -4, -2, 11, 23.741881),
]


@pytest.mark.parametrize(
    ("integrand", "start", "end", "sample_count", "result"), INTEGRALS
)
def test_integral(integrand, start, end, sample_count, result):
    x = samples(start, end, sample_count)
    f = np.asarray(list(map(integrand, x)))
    i = integral(f, x)
    assert i == pytest.approx(result)


@pytest.mark.parametrize(
    ("integrand", "start", "end", "sample_count", "result"), INTEGRALS
)
def test_integral_equals_np_trapz(integrand, start, end, sample_count, result):
    x = samples(start, end, sample_count)
    f = np.asarray(list(map(integrand, x)))
    i = integral(f, x)
    assert i == pytest.approx(np.trapz(f, x))


LINE_INTEGRALS = [
    (lambda t: [t, 0], lambda t: [0.5 * t, 0], 0, 1, 10, 0.25),
    (lambda t: [math.sin(t), 0], lambda t: [t, t], 1, 3, 50, 1.5300823),
    (lambda t: [t**2, -t], lambda t: [math.cos(t), 0], -4, 7, 500, 29.33303),
]


@pytest.mark.parametrize(
    ("integrand", "path", "start", "end", "sample_count", "result"), LINE_INTEGRALS
)
def test_line_integral(integrand, path, start, end, sample_count, result):
    t = samples(start, end, sample_count)
    x = np.asarray(list(map(path, t)))
    f = np.asarray(list(map(integrand, t)))
    integral = vector_line_integral(f, x)
    assert integral == pytest.approx(result)


@pytest.mark.parametrize(
    ("integrand", "path", "start", "end", "sample_count", "result"), LINE_INTEGRALS
)
def test_cum_line_integral(integrand, path, start, end, sample_count, result):
    t = samples(start, end, sample_count)
    x = np.asarray(list(map(path, t)))
    f = np.asarray(list(map(integrand, t)))
    integral = vector_line_integral(f, x)
    cum_integral = cumulative_vector_line_integral(f, x)
    assert cum_integral[..., -1] == pytest.approx(integral)
