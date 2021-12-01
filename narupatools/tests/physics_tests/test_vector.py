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
    outer_product,
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


def test_zero_vector():
    assert zero_vector() == pytest.approx([0.0, 0.0, 0.0])


def test_magnitude_zero():
    assert magnitude(zero_vector()) == pytest.approx(0.0)


def test_sqr_magnitude_zero():
    assert sqr_magnitude(zero_vector()) == pytest.approx(0.0)


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


def outer_calc(vec1, vec2):
    vec1 = np.asarray(vec1)
    vec2 = np.asarray(vec2)
    if len(vec1.shape) != 1 or vec1.shape != vec2.shape:
        raise ValueError(f"Cannot take outer product of {vec1} and {vec2}")
    out = np.zeros((len(vec1), len(vec2)))
    for i in range(len(vec1)):
        for j in range(len(vec2)):
            out[i, j] = vec1[i] * vec2[j]
    return pytest.approx(out)


def dot_calc(vec1, vec2):
    vec1 = np.asarray(vec1)
    vec2 = np.asarray(vec2)
    if len(vec1.shape) != 1 or vec1.shape != vec2.shape:
        raise ValueError(f"Cannot take dot product of {vec1} and {vec2}")
    dot = 0.0
    for a, b in zip(vec1, vec2):
        dot += a * b
    return dot


def cross_calc(vec1, vec2):
    vec1 = np.asarray(vec1)
    vec2 = np.asarray(vec2)
    if len(vec1.shape) != 1 or vec1.shape != vec2.shape or vec1.shape != (3,):
        raise ValueError(f"Cannot take cross product of {vec1} and {vec2}")
    return pytest.approx(
        np.array(
            [
                vec1[1] * vec2[2] - vec1[2] * vec2[1],
                vec1[2] * vec2[0] - vec1[0] * vec2[2],
                vec1[0] * vec2[1] - vec1[1] * vec2[0],
            ]
        )
    )


def sqr_distance_calc(vec1, vec2):
    vec1 = np.asarray(vec1)
    vec2 = np.asarray(vec2)
    if len(vec1.shape) != 1 or vec1.shape != vec2.shape:
        raise ValueError(f"Cannot take distance of {vec1} and {vec2}")
    dist2 = 0.0
    for a, b in zip(vec1, vec2):
        dist2 += (a - b) ** 2
    return pytest.approx(dist2)


def distance_calc(vec1, vec2):
    vec1 = np.asarray(vec1)
    vec2 = np.asarray(vec2)
    if len(vec1.shape) != 1 or vec1.shape != vec2.shape:
        raise ValueError(f"Cannot take distance of {vec1} and {vec2}")
    dist2 = 0.0
    for a, b in zip(vec1, vec2):
        dist2 += (a - b) ** 2
    return pytest.approx(math.sqrt(dist2))


def angle_calc(vec1, vec2):
    vec1 = np.asarray(vec1)
    vec2 = np.asarray(vec2)
    if len(vec1.shape) != 1 or vec1.shape != vec2.shape:
        raise ValueError(f"Cannot take distance of {vec1} and {vec2}")
    mag1 = 0.0
    dot = 0.0
    mag2 = 0.0
    for a, b in zip(vec1, vec2):
        mag1 += a * a
        mag2 += b * b
        dot += a * b
    mag1 = math.sqrt(mag1)
    mag2 = math.sqrt(mag2)
    return pytest.approx(math.acos(dot / (mag1 * mag2)))


def scalar_projection_calc(vec1, vec2):
    vec1 = np.asfarray(vec1)
    vec2 = np.asfarray(vec2)
    dot = 0.0
    mag2 = 0.0
    for a, b in zip(vec1, vec2):
        dot += a * b
        mag2 += b * b
    return dot / math.sqrt(mag2)


def vector_projection_calc(vec1, vec2):
    vec1 = np.asfarray(vec1)
    vec2 = np.asfarray(vec2)
    dot = 0.0
    mag2 = 0.0
    for a, b in zip(vec1, vec2):
        dot += a * b
        mag2 += b * b
    return pytest.approx(dot / mag2 * vec2)


def vector_rejection_calc(vec1, vec2):
    vec1 = np.asfarray(vec1)
    vec2 = np.asfarray(vec2)
    dot = 0.0
    mag2 = 0.0
    for a, b in zip(vec1, vec2):
        dot += a * b
        mag2 += b * b
    return pytest.approx(vec1 - dot / mag2 * vec2)


@pytest.mark.parametrize(
    ("op", "op_test"),
    [
        (outer_product, outer_calc),
        (dot_product, dot_calc),
        (cross_product, cross_calc),
        (sqr_distance, sqr_distance_calc),
        (distance, distance_calc),
        (angle, angle_calc),
        (vector_projection, vector_projection_calc),
        (vector_rejection, vector_rejection_calc),
    ],
)
def test_vector_op(op, op_test):
    vec1 = np.random.random((3,))
    vec2 = np.random.random((3,))
    assert op(vec1, vec2) == op_test(vec1, vec2)


@pytest.mark.parametrize(
    ("op", "op_test", "shape"),
    [
        (outer_product, outer_calc, (3, 3)),
        (dot_product, dot_calc, ()),
        (cross_product, cross_calc, (3,)),
        (sqr_distance, sqr_distance_calc, ()),
        (distance, distance_calc, ()),
        (angle, angle_calc, ()),
        (vector_projection, vector_projection_calc, (3,)),
        (vector_rejection, vector_rejection_calc, (3,)),
    ],
)
def test_vector_op_left_array(op, op_test, shape):
    vecs1 = np.random.random((10, 3))
    vec2 = np.random.random((3,))
    result = op(vecs1, vec2)
    assert result.shape == (10, *shape)
    for i in range(10):
        assert result[i] == op_test(vecs1[i], vec2)


@pytest.mark.parametrize(
    ("op", "op_test", "shape"),
    [
        (outer_product, outer_calc, (3, 3)),
        (dot_product, dot_calc, ()),
        (cross_product, cross_calc, (3,)),
        (sqr_distance, sqr_distance_calc, ()),
        (distance, distance_calc, ()),
        (angle, angle_calc, ()),
        (vector_projection, vector_projection_calc, (3,)),
        (vector_rejection, vector_rejection_calc, (3,)),
    ],
)
def test_vector_op_right_array(op, op_test, shape):
    vec1 = np.random.random((3))
    vecs2 = np.random.random((10, 3))
    result = op(vec1, vecs2)
    assert result.shape == (10, *shape)
    for i in range(10):
        assert result[i] == op_test(vec1, vecs2[i])


@pytest.mark.parametrize(
    ("op", "op_test", "shape"),
    [
        (outer_product, outer_calc, (3, 3)),
        (dot_product, dot_calc, ()),
        (cross_product, cross_calc, (3,)),
        (sqr_distance, sqr_distance_calc, ()),
        (distance, distance_calc, ()),
        (angle, angle_calc, ()),
        (vector_projection, vector_projection_calc, (3,)),
        (vector_rejection, vector_rejection_calc, (3,)),
    ],
)
def test_vector_op_both_array(op, op_test, shape):
    vecs1 = np.random.random((10, 3))
    vecs2 = np.random.random((10, 3))
    result = op(vecs1, vecs2)
    assert result.shape == (10, *shape)
    for i in range(10):
        assert result[i] == op_test(vecs1[i], vecs2[i])


def normalized_calc(vec):
    return vec / magnitude_calc(vec)


def sqr_magnitude_calc(vec):
    mag2 = 0.0
    for i in vec:
        mag2 += i * i
    return mag2


def magnitude_calc(vec):
    return math.sqrt(sqr_magnitude_calc(vec))


@pytest.mark.parametrize(
    ("op", "op_test"),
    [
        (normalized, normalized_calc),
        (magnitude, magnitude_calc),
        (sqr_magnitude, sqr_magnitude_calc),
    ],
)
def test_vector_1op(op, op_test):
    vec1 = np.random.random((3,))
    assert op(vec1) == pytest.approx(op_test(vec1))


@pytest.mark.parametrize(
    ("op", "op_test", "shape"),
    [
        (normalized, normalized_calc, (3,)),
        (magnitude, magnitude_calc, ()),
        (sqr_magnitude, sqr_magnitude_calc, ()),
    ],
)
def test_vector_1op_array(op, op_test, shape):
    vecs = np.random.random((10, 3))
    result = op(vecs)
    assert result.shape == (10, *shape)
    for i in range(10):
        assert result[i] == pytest.approx(op_test(vecs[i]))
