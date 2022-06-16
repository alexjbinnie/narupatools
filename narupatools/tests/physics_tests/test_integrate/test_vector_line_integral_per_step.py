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

import numpy as np

from narupatools.physics.integrate import vector_line_integral_per_step


def test_shape():
    samples = 30
    components = 3
    integrand = np.random.random((samples, components))
    x = np.random.random((samples, components))
    integral = vector_line_integral_per_step(integrand, x)
    assert integral.shape == (samples - 1,)


def test_axis_shape():
    samples = 30
    grouped = 7
    components = 3
    integrand = np.random.random((samples, grouped, components))
    x = np.random.random((samples, grouped, components))
    integral = vector_line_integral_per_step(integrand, x, axis=-2)
    assert integral.shape == (samples - 1, grouped)


def test_variable_group_shape():
    samples = 30
    grouped = 7
    components = 3
    integrand = np.random.random((samples, components))
    x = np.random.random((grouped, samples, components))
    integral = vector_line_integral_per_step(integrand, x)
    assert integral.shape == (grouped, samples - 1)


def test_integrand_group_shape():
    samples = 30
    grouped = 7
    components = 3
    integrand = np.random.random((grouped, samples, components))
    x = np.random.random((samples, components))
    integral = vector_line_integral_per_step(integrand, x)
    assert integral.shape == (grouped, samples - 1)


def test_both_group_shape():
    samples = 30
    grouped = 7
    components = 3
    integrand = np.random.random((grouped, samples, components))
    x = np.random.random((grouped, samples, components))
    integral = vector_line_integral_per_step(integrand, x)
    assert integral.shape == (grouped, samples - 1)
