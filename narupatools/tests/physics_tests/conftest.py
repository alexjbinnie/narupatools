import random

import pytest


@pytest.fixture(params=range(20))
def seed(request):
    random.seed(request.param)
    return request.param
