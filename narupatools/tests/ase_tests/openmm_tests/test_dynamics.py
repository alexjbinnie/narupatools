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
from test_classes.dynamics import VillinDynamicsTests

from narupatools.ase import ASEDynamics
from narupatools.ase.openmm import openmm_simulation_to_ase_molecular_dynamics


@pytest.fixture
def villin_molecular_dynamics(villin_openmm_simulation):
    return openmm_simulation_to_ase_molecular_dynamics(villin_openmm_simulation)


@pytest.fixture
def villin_dynamics(villin_molecular_dynamics):
    return ASEDynamics(villin_molecular_dynamics)


@pytest.mark.dynamics
class TestVillinASEOpenMMDynamics(VillinDynamicsTests):
    @pytest.fixture(autouse=True)
    def dynamics(self, villin_dynamics):
        return villin_dynamics

