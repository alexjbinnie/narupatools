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
from ase.md import VelocityVerlet
from simtk.openmm import System
from simtk.openmm.app import Element, Simulation, Topology
from test_classes.single_carbon import SingleCarbonSystemTests

from narupatools.ase import ASEDynamics, UnitsASE
from narupatools.ase.openmm import openmm_simulation_to_ase_atoms
from narupatools.core import UnitsNarupa
from narupatools.openmm.dynamics import OpenMMDynamics
from narupatools.openmm.integrators import velocity_verlet_integrator
from narupatools.physics.vector import vector

_NarupaToASE = UnitsNarupa >> UnitsASE


@pytest.fixture
def simulation():
    system = System()
    system.addParticle(12.000)
    topology = Topology()
    chain = topology.addChain(0)
    residue = topology.addResidue("RES", chain)
    topology.addAtom("C", Element.getBySymbol("C"), residue)
    integrator = velocity_verlet_integrator(0.01)
    simulation = Simulation(topology, system, integrator)
    simulation.context.setPositions([vector(5.0, 5.0, 5.0)])
    return simulation


class TestOpenMMSingleCarbonSystem(SingleCarbonSystemTests):
    @pytest.fixture
    def dynamics(self, simulation):
        return OpenMMDynamics(simulation)


class TestASEOpenMMSingleCarbonSystem(SingleCarbonSystemTests):
    @pytest.fixture
    def dynamics(self, simulation):
        atoms = openmm_simulation_to_ase_atoms(simulation)
        verlet = VelocityVerlet(atoms, timestep=0.01 * _NarupaToASE.time)
        dynamics = ASEDynamics(verlet)
        return dynamics
