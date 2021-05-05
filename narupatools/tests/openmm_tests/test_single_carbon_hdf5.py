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

import pytest
from ase.md import Langevin
from simtk.openmm import LangevinIntegrator, System
from simtk.openmm.app import Element, Simulation, Topology

from narupatools.ase import ASEDynamics, UnitsASE
from narupatools.ase.openmm import openmm_simulation_to_ase_atoms
from narupatools.core import UnitsNarupa
from narupatools.openmm.dynamics import OpenMMDynamics
from narupatools.physics.vector import vector
from test_classes.single_carbon_hdf5 import SingleCarbonHDF5Tests

_NarupaToASE = UnitsNarupa >> UnitsASE


@pytest.fixture
def simulation():
    system = System()
    system.addParticle(12.000)
    topology = Topology()
    chain = topology.addChain(0)
    residue = topology.addResidue("RES", chain)
    topology.addAtom("C", Element.getBySymbol("C"), residue)
    integrator = LangevinIntegrator(300, 0.01, 0.01)
    simulation = Simulation(topology, system, integrator)
    simulation.context.setPositions([vector(5.0, 5.0, 5.0)])
    return simulation


class TestOpenMMSingleCarbonHDF5(SingleCarbonHDF5Tests):
    @pytest.fixture
    def dynamics(self, simulation):
        return OpenMMDynamics(simulation)


class TestASEOpenMMSingleCarbonHDF5(SingleCarbonHDF5Tests):
    @pytest.fixture
    def dynamics(self, simulation):
        atoms = openmm_simulation_to_ase_atoms(simulation)
        langevin = Langevin(
            atoms,
            timestep=0.01 * _NarupaToASE.time,
            temperature_K=300,
            friction=0.01 / _NarupaToASE.time,
        )
        dynamics = ASEDynamics(langevin)
        return dynamics
