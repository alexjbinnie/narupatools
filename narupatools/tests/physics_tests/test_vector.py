import math

import numpy as np
import pytest

from narupatools.physics._quaternion import quaternion
from narupatools.physics.random import (
    random_float,
    random_integer,
    random_quaternion,
    random_vector,
)
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
    return random_float(minimum=-100.0, maximum=100.0)


@pytest.fixture
def y():
    return random_float(minimum=-100.0, maximum=100.0)


@pytest.fixture
def z():
    return random_float(minimum=-100.0, maximum=100.0)


@pytest.fixture
def vec3(seed):
    return random_vector(max_magnitude=100.0)


@pytest.fixture
def vec3_list(seed):
    return [
        random_vector(max_magnitude=100.0)
        for _ in range(random_integer(minimum=5, maximum=50))
    ]


@pytest.fixture
def vec3_array(vec3_list):
    return np.array(vec3_list)


@pytest.fixture
def quat(seed):
    return random_quaternion()


@pytest.fixture
def quat_list(seed):
    return [random_quaternion() for _ in range(random_integer(minimum=5, maximum=50))]


@pytest.fixture
def quat_array(quat_list):
    return np.array(quat_list, dtype=quaternion)


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


def test_magnitude(vec3):
    assert magnitude(vec3) == pytest.approx(np.linalg.norm(vec3))


def test_magnitude_zero():
    assert magnitude(zero_vector()) == pytest.approx(0.0)


def test_sqr_magnitude(vec3):
    assert sqr_magnitude(vec3) == pytest.approx(np.dot(vec3, vec3))


def test_sqr_magnitude_zero():
    assert sqr_magnitude(zero_vector()) == pytest.approx(0.0)


def test_normalize_vec(vec3):
    n = normalized(vec3)
    assert n[0] * n[0] + n[1] * n[1] + n[2] * n[2] == pytest.approx(1.0)


def test_normalize_vec3_array(vec3_array):
    vec3_array = normalized(vec3_array)
    for n in vec3_array:
        assert n[0] * n[0] + n[1] * n[1] + n[2] * n[2] == pytest.approx(1.0)


def test_normalize_vec3_list(vec3_list):
    vec3_list = normalized(vec3_list)
    for n in vec3_list:
        assert n[0] * n[0] + n[1] * n[1] + n[2] * n[2] == pytest.approx(1.0)


def test_normalize_quaternion(quat):
    n = normalized(quat)
    assert n.w * n.w + n.x * n.x + n.y * n.y + n.z * n.z == pytest.approx(1.0)


def test_normalize_quat_array(quat_array):
    quat_array = normalized(quat_array)
    for n in quat_array:
        assert n.w * n.w + n.x * n.x + n.y * n.y + n.z * n.z == pytest.approx(1.0)


def test_normalize_quat_list(quat_list):
    quat_list = normalized(quat_list)
    for n in quat_list:
        assert n.w * n.w + n.x * n.x + n.y * n.y + n.z * n.z == pytest.approx(1.0)


def test_vector_normalized_zero_vector():
    assert normalized(zero_vector()) == pytest.approx(zero_vector())


def test_cross_product(vec1, vec2):
    assert cross_product(vec1, vec2) == pytest.approx(np.cross(vec1, vec2))


def test_cross_product_matrix(vec1, vec2):
    mat = cross_product_matrix(vec1)
    assert mat.shape == (3, 3)
    assert mat @ vec2 == pytest.approx(np.cross(vec1, vec2))


def test_right_cross_product_matrix(vec1, vec2):
    mat = right_cross_product_matrix(vec2)
    assert mat.shape == (3, 3)
    assert mat @ vec1 == pytest.approx(np.cross(vec1, vec2))


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
    assert distance(vec1, vec2) == pytest.approx(magnitude(offset))


def test_sqr_distance(vec1, vec2):
    offset = np.subtract(vec1, vec2)
    assert sqr_distance(vec1, vec2) == pytest.approx(sqr_magnitude(offset))


def test_angle_zero(vec3):
    with pytest.raises(ValueError):  # noqa: PT011
        _ = angle(vec3, zero_vector())
    with pytest.raises(ValueError):  # noqa: PT011
        _ = angle(zero_vector(), vec3)
    with pytest.raises(ValueError):  # noqa: PT011
        _ = angle(zero_vector(), zero_vector())


def test_projection_magnitude(vec1, vec2):
    projection = vector_projection(vec1, vec2)
    assert magnitude(projection) == pytest.approx(
        abs(magnitude(vec1) * math.cos(angle(vec1, vec2)))
    )
