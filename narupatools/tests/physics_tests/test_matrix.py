import numpy as np
import pytest

from narupatools.physics.matrix import kronecker_delta, zero_matrix


def test_zero_matrix():
    assert zero_matrix() == pytest.approx(
        np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    )


@pytest.mark.parametrize("i", [0, 1, 3])
@pytest.mark.parametrize("j", [0, 1, 3])
def test_kronecker_delta(i, j):
    assert kronecker_delta(i, j) == (1 if i == j else 0)
