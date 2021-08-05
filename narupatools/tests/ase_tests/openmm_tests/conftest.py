import pytest
from testing import add_mark


def pytest_collection_modifyitems(items):
    add_mark(filename=__file__, mark=pytest.mark.openmm, items=items)
