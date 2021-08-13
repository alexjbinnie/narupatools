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

from io import BytesIO

import pytest
from simtk.openmm import LangevinIntegrator
from simtk.openmm.app import ForceField, HBonds, PDBFile, Simulation
from simtk.openmm.app.forcefield import PME
from simtk.unit import kelvin, nanometer, picosecond, picoseconds
from test_classes.dynamics import VillinDynamicsTests

from narupatools.openmm.dynamics import OpenMMDynamics


@pytest.fixture(scope="module")
def villin_pdbfile(villin_pdb_filename) -> PDBFile:
    return PDBFile(villin_pdb_filename)


@pytest.fixture(scope="module")
def villin_simulation_original(villin_pdbfile) -> Simulation:
    forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
    system = forcefield.createSystem(
        villin_pdbfile.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1 * nanometer,
        constraints=HBonds,
    )
    integrator = LangevinIntegrator(300 * kelvin, 1 / picosecond, 0.005 * picoseconds)
    simulation = Simulation(villin_pdbfile.topology, system, integrator)
    simulation.context.setPositions(villin_pdbfile.positions)
    return simulation


@pytest.fixture(scope="module")
def villin_simulation_checkpoint(villin_simulation_original) -> bytes:
    with BytesIO() as bytesio:
        villin_simulation_original.saveCheckpoint(bytesio)
        return bytesio.getvalue()


@pytest.fixture
def villin_simulation(villin_simulation_original, villin_simulation_checkpoint):
    with BytesIO(villin_simulation_checkpoint) as bytesio:
        villin_simulation_original.loadCheckpoint(bytesio)
        return villin_simulation_original


@pytest.fixture
def villin_dynamics(villin_simulation):
    return OpenMMDynamics(villin_simulation)


@pytest.mark.dynamics
class TestVillinOpenMMDynamics(VillinDynamicsTests):
    @pytest.fixture(autouse=True)
    def dynamics(self, villin_dynamics):
        return villin_dynamics

    def test_reset_works(self, dynamics):
        positions = dynamics.positions
        velocities = dynamics.velocities
        forces = dynamics.forces
        dynamics.run(100)
        dynamics.reset()
        assert dynamics.positions == pytest.approx(positions, rel=1e-3)
        assert dynamics.velocities == pytest.approx(velocities, rel=1e-3)
        assert dynamics.forces == pytest.approx(forces, rel=1e-1)
