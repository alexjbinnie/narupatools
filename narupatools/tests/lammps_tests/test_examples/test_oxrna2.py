import pytest

pytest.importorskip("lammps")


from narupatools.lammps import LAMMPSSimulation


def test_simulation():
    _ = LAMMPSSimulation.from_file("in.duplex4")
