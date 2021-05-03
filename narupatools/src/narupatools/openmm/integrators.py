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
#
# Some code taken from the documentation of OpenMM.
# Copyright (c) 2011-2020 Stanford University and the Authors.

"""Custom integrators for use with OpenMM."""

from simtk.openmm import CustomIntegrator


def velocity_verlet_integrator(timestep: float) -> CustomIntegrator:
    """
    Velocity-verlet integrator, using RATTLE to implement constraints.

    This is taken from the documentation for the OpenMM CustomIntegrator class.

    :param timestep: Timestep in picoseconds.
    :return: Custom integrator that implements velocity-verlet.
    """
    integrator = CustomIntegrator(timestep)
    integrator.addPerDofVariable("x1", 0)
    integrator.addUpdateContextState()
    integrator.addComputePerDof("v", "v+0.5*dt*f/m")
    integrator.addComputePerDof("x", "x+dt*v")
    integrator.addComputePerDof("x1", "x")
    integrator.addConstrainPositions()
    integrator.addComputePerDof("v", "v+0.5*dt*f/m+(x-x1)/dt")
    integrator.addConstrainVelocities()
    return integrator
