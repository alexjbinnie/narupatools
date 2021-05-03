# This file is part of narupatools (https://github.com/alexjbinnie/narupatools).
# Copyright (c) Alex Jamieson-Binnie. All rights reserved.
#
# narupatools is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# narupatools is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with narupatools.  If not, see <http://www.gnu.org/licenses/>.

import os
import sys
from io import BytesIO
from pathlib import Path

import ase
import pytest
from MDAnalysis import Universe
from ase import Atoms
from ase.md import Langevin
from ase.units import second
from simtk.openmm import LangevinIntegrator
from simtk.openmm.app import ForceField, HBonds, PDBFile, PME, Simulation
from simtk.unit import kelvin, nanometer, picosecond, picoseconds

from narupatools.ase import ASEDynamics
from narupatools.ase import NullCalculator

sys.path.append(os.path.join(os.path.dirname(__file__), "testing"))


@pytest.fixture(scope="session")
def neuraminidase_pdb_filename():
    return str(Path(__file__).parent / "files/2qwk.pdb")


@pytest.fixture(scope="session")
def ethane_sdf_filename():
    return str(Path(__file__).parent / "files/ethane.sdf")


@pytest.fixture(scope="session")
def villin_pdb_filename():
    return str(Path(__file__).parent / "files/villin.pdb")


@pytest.fixture(scope="session")
def nanotube_xml_filename():
    return str(Path(__file__).parent / "files/nanotube.xml")


@pytest.fixture(scope="session")
def villin_openmm_pdbfile(villin_pdb_filename) -> PDBFile:
    return PDBFile(villin_pdb_filename)


@pytest.fixture(scope="session")
def neuraminidase_openmm_xml_filename() -> str:
    return str(Path(__file__).parent / "files/neuraminidase.xml")


@pytest.fixture(scope="session")
def villin_openmm_topology(villin_pdbfile):
    return villin_pdbfile.topology


@pytest.fixture(scope="session")
def villin_openmm_system(villin_openmm_pdbfile):
    forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
    return forcefield.createSystem(
        villin_openmm_pdbfile.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1 * nanometer,
        constraints=HBonds,
    )


@pytest.fixture(scope="session")
def villin_openmm_simulation_original(villin_openmm_pdbfile) -> Simulation:
    forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
    system = forcefield.createSystem(
        villin_openmm_pdbfile.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1 * nanometer,
        constraints=HBonds,
    )
    integrator = LangevinIntegrator(300 * kelvin, 1 / picosecond, 0.005 * picoseconds)
    simulation = Simulation(villin_openmm_pdbfile.topology, system, integrator)
    simulation.context.setPositions(villin_openmm_pdbfile.positions)
    return simulation


@pytest.fixture(scope="session")
def villin_openmm_simulation_checkpoint(villin_openmm_simulation_original) -> bytes:
    with BytesIO() as bytesio:
        villin_openmm_simulation_original.saveCheckpoint(bytesio)
        return bytesio.getvalue()


@pytest.fixture
def villin_openmm_simulation(
    villin_openmm_simulation_original, villin_openmm_simulation_checkpoint
):
    with BytesIO(villin_openmm_simulation_checkpoint) as bytesio:
        villin_openmm_simulation_original.loadCheckpoint(bytesio)
        return villin_openmm_simulation_original


@pytest.fixture(scope="session")
def villin_ase_atoms_readonly(villin_pdb_filename) -> Atoms:
    return ase.io.read(villin_pdb_filename)  # type: ignore


@pytest.fixture
def villin_ase_atoms(villin_ase_atoms_readonly) -> Atoms:
    return villin_ase_atoms_readonly.copy()  # type: ignore


@pytest.fixture
def villin_ase_langevin_dynamics(villin_ase_atoms) -> Langevin:
    atoms = villin_ase_atoms.copy()
    atoms.set_calculator(NullCalculator())
    dynamics = Langevin(
        atoms,
        timestep=0.005 * (1e-12 * second),
        friction=1.0 / (1e-12 * second),
        temperature_K=300,
    )
    return dynamics


@pytest.fixture
def villin_ase_simulation_dynamics(villin_ase_langevin_dynamics) -> ASEDynamics:
    return ASEDynamics(villin_ase_langevin_dynamics)


@pytest.fixture(scope="session")
def villin_mda_universe_readonly(villin_pdb_filename) -> Universe:
    return Universe(villin_pdb_filename, guess_bonds=True)


pytest.register_assert_rewrite("tests.base_tests.converter")
pytest.register_assert_rewrite("tests.base_tests.dynamics")


# Skip OpenMM exceptions, which occur on CI but not locally.
def pytest_runtest_makereport(item, call):
    from _pytest.runner import pytest_runtest_makereport as _pytest_runtest_makereport
    from simtk.openmm import OpenMMException

    tr = _pytest_runtest_makereport(item, call)

    if call.excinfo is not None and call.excinfo.type == OpenMMException:
        tr.outcome = "skipped"
        tr.wasxfail = "reason: OpenMMException occurs on CLI"

    return tr
