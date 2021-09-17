import numpy as np
import pytest
from ase import Atoms
from ase.io import read
from MDAnalysis import Universe
from narupa.trajectory import FrameData

from narupatools.ase import NullCalculator
from narupatools.frame import convert


@pytest.fixture(scope="module")
def ase_atoms(neuraminidase_pdb_filename):
    atoms = read(neuraminidase_pdb_filename)
    atoms.set_calculator(NullCalculator())
    return atoms


@pytest.fixture(scope="module")
def mda_universe(neuraminidase_pdb_filename):
    return Universe(neuraminidase_pdb_filename, guess_bonds=True)


def test_empty_frame():
    _ = convert(FrameData(), Atoms)
    _ = convert(FrameData(), Universe)


def test_to_from_frame(ase_atoms):
    frame = convert(ase_atoms, FrameData)
    ase_atoms2 = convert(frame, Atoms)
    assert ase_atoms.get_positions() == pytest.approx(ase_atoms2.get_positions())
    assert np.array_equal(ase_atoms.symbols, ase_atoms2.symbols)
    assert ase_atoms.get_momenta() == pytest.approx(ase_atoms2.get_momenta())


def test_to_from_universe(mda_universe):
    frame = convert(mda_universe, FrameData)
    mda_universe2 = convert(frame, Universe)
    assert mda_universe.atoms.positions == pytest.approx(mda_universe2.atoms.positions)


def test_ase_atoms_to_mda_universe(ase_atoms, mda_universe):
    mda_universe2 = convert(ase_atoms, Universe)
    assert mda_universe.atoms.positions == pytest.approx(mda_universe2.atoms.positions)
