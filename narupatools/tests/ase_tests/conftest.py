import numpy as np
import pytest
from ase import Atoms
from ase.io import read


@pytest.fixture(scope="session")
def ethane_atoms_readonly(ethane_sdf_filename) -> Atoms:
    return read(ethane_sdf_filename)  # type: ignore[return-value]


@pytest.fixture
def ethane_atoms(ethane_atoms_readonly) -> Atoms:
    return ethane_atoms_readonly.copy()  # type: ignore[no-any-return]


@pytest.fixture
def single_carbon_atoms() -> Atoms:
    return Atoms(
        symbols=["C"], masses=np.array([12.0]), positions=np.array([[0.0, 0.0, 0.0]])
    )
