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

"""Classes and methods for interfacing with the OpenMM package."""

from ._converter import openmm_state_to_frame, openmm_topology_to_frame
from ._dynamics import OpenMMDynamics
from ._integrators import VelocityVerletIntegrator
from ._serializer import deserialize_simulation, serialize_simulation
from ._subset import simulation_subset
from ._units import UnitsOpenMM
from ._simulation import OpenMMSimulation

__all__ = [
    "openmm_state_to_frame",
    "openmm_topology_to_frame",
    "UnitsOpenMM",
    "serialize_simulation",
    "deserialize_simulation",
    "OpenMMDynamics",
    "VelocityVerletIntegrator",
    "simulation_subset",
    "OpenMMSimulation",
]
