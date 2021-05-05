import math

import numpy as np
import pytest

from narupatools.physics.random import random_scalar, random_vector
from narupatools.physics.vector import (
    angle,
    cross_product,
    cross_product_matrix,
    distance,
    dot_product,
    magnitude,
    normalized,
    right_cross_product_matrix,
    sqr_distance,
    sqr_magnitude,
    vector,
    vector_projection,
    vector_rejection,
    zero_vector,
)


@pytest.fixture
def x():
    return random_scalar(min=-100.0, max=100.0)


@pytest.fixture
def y():
    return random_scalar(min=-100.0, max=100.0)


@pytest.fixture
def z():
    return random_scalar(min=-100.0, max=100.0)


@pytest.fixture
def vec():
    return random_vector(max_magnitude=100.0)


@pytest.fixture
def vec1():
    return random_vector(max_magnitude=100.0)


@pytest.fixture
def vec2():
    return random_vector(max_magnitude=100.0)


def test_vector(x, y, z):
    vec = vector(x, y, z)
    assert vec == pytest.approx(np.array([x, y, z]))
    assert isinstance(vec, np.ndarray)
    assert vec.dtype == float


def test_dot_product(vec1, vec2):
    assert dot_product(vec1, vec2) == pytest.approx(np.dot(vec1, vec2))


def test_zero_vector():
    assert zero_vector() == pytest.approx([0.0, 0.0, 0.0])


def test_magnitude(vec):
    assert magnitude(vec) == pytest.approx(np.linalg.norm(vec))


def test_magnitude_zero():
    assert magnitude(zero_vector()) == pytest.approx(0.0)


def test_sqr_magnitude(vec):
    assert sqr_magnitude(vec) == pytest.approx(np.dot(vec, vec))


def test_sqr_magnitude_zero():
    assert sqr_magnitude(zero_vector()) == pytest.approx(0.0)


def test_vector_normalized(vec):
    n = normalized(vec)
    assert magnitude(n) == pytest.approx(1.0)
    assert dot_product(n, vec) == pytest.approx(magnitude(vec))


def test_vector_normalized_zero_vector():
    assert normalized(zero_vector()) == pytest.approx(zero_vector())


def test_cross_product(vec1, vec2):
    assert cross_product(vec1, vec2) == pytest.approx(np.cross(vec1, vec2))


def test_cross_product_matrix(vec1, vec2):
    mat = cross_product_matrix(vec1)
    assert mat.shape == (3, 3)
    assert np.matmul(mat, vec2) == pytest.approx(np.cross(vec1, vec2))


def test_right_cross_product_matrix(vec1, vec2):
    mat = right_cross_product_matrix(vec2)
    assert mat.shape == (3, 3)
    assert np.matmul(mat, vec1) == pytest.approx(np.cross(vec1, vec2))


def test_projection_rejection_sum(vec1, vec2):
    projection = vector_projection(vec1, vec2)
    rejection = vector_rejection(vec1, vec2)
    assert vec1 == pytest.approx(projection + rejection)


def test_projection_rejection_perpindicular(vec1, vec2):
    projection = vector_projection(vec1, vec2)
    rejection = vector_rejection(vec1, vec2)
    assert dot_product(projection, rejection) == pytest.approx(0.0)


def test_angle(vec1, vec2):
    norm1 = magnitude(vec1)
    norm2 = magnitude(vec2)
    a = angle(vec1, vec2)
    assert math.cos(a) * norm1 * norm2 == pytest.approx(dot_product(vec1, vec2))


def test_distance(vec1, vec2):
    offset = np.subtract(vec1, vec2)
    assert distance(vec1, vec2) == magnitude(offset)


def test_sqr_distance(vec1, vec2):
    offset = np.subtract(vec1, vec2)
    assert sqr_distance(vec1, vec2) == sqr_magnitude(offset)


def test_angle_zero(vec):
    with pytest.raises(ValueError):  # noqa: PT011
        _ = angle(vec, zero_vector())
    with pytest.raises(ValueError):  # noqa: PT011
        _ = angle(zero_vector(), vec)
    with pytest.raises(ValueError):  # noqa: PT011
        _ = angle(zero_vector(), zero_vector())


def test_projection_magnitude(vec1, vec2):
    projection = vector_projection(vec1, vec2)
    assert magnitude(projection) == pytest.approx(
        abs(magnitude(vec1) * math.cos(angle(vec1, vec2)))
    )
