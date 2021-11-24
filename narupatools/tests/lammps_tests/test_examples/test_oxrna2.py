import pytest

pytest.importorskip("lammps")


from narupatools.lammps import LAMMPSSimulation


def test_simulation():
    simulation = LAMMPSSimulation.from_file("in.duplex4")
    simulation.run(100)
