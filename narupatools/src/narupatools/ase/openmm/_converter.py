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
import functools
import logging
from collections import Callable

from ase.atoms import Atoms
from ase.md import Langevin, VelocityVerlet
from ase.md.md import MolecularDynamics
from simtk.openmm import LangevinIntegrator, VerletIntegrator
from simtk.openmm.app import Simulation
from simtk.unit import picoseconds

import narupatools.ase.openmm._calculator as omm_calculator
from narupatools.ase import UnitsASE
from narupatools.frame import convert
from narupatools.openmm import UnitsOpenMM
from narupatools.physics.units import UnitsNarupa, pico, second

_OpenMMToASE = UnitsOpenMM >> UnitsASE
_NarupaToASE = UnitsNarupa >> UnitsASE
_ASEToNarupa = UnitsASE >> UnitsNarupa

DEFAULT_LANGEVIN_FRICTION = 10.0 / (pico * second)  # Friction in per picosecond

INTEGRATOR_VERLET_MESSAGE = (
    "Running OpenMM simulation that uses OpenMM's leap-frog Verlet integrator. ASE will use "
    "a velocity Verlet integrator with the same timestep."
)
INTEGRATOR_NOT_LANGEVIN_MESSAGE = (
    "Running OpenMM simulation that was not setup with Langevin Integrator. A Langevin "
    "integrator with friction 0.01 fs^{-1} will be used."
)
CONSTRAINTS_UNSUPPORTED_MESSAGE = (
    "The simulation contains constraints which will be ignored by this runner!"
)


def openmm_simulation_to_ase_atoms(simulation: Simulation, /) -> Atoms:
    """
    Generate an ASE atoms representation of the OpenMM simulation.

    The :class:`~ase.atoms.Atoms` object will have a
    :class:`~narupatools.ase.openmm.OpenMMCalculator` which calculates the forces and
    energy based on the simulation.

    :param simulation: OpenMM simulation to convert to :class:`~ase.atoms.Atoms`.
    :return: ASE atoms object, with positions and chemical symbols read from the current
             state of the OpenMM simulation.
    """
    atoms = convert(simulation, Atoms)
    atoms.set_pbc(simulation.system.usesPeriodicBoundaryConditions())

    atoms.calc = omm_calculator.OpenMMCalculator(simulation)

    return atoms


def openmm_simulation_to_ase_molecular_dynamics(
    simulation: Simulation,
) -> MolecularDynamics:
    """
    Convert an OpenMM simulation to an ASE simulation.

    This both converts the system itself to an ASE atoms object
    and converts the simulation to an ASE Langevin integrator. The timestep and
    temperature are read from the OpenMM simulation's integrator. If the temperature is
    not present, a default of 300 kelvin is used. If the OpenMM integrator is Langevin,
    then  the friction is also copied. If no friction is present, a default friction of
    10 ps^{-1} is used.

    The ASE atoms object has a calculator which will call the original underlying
    simulation to calculate forces.

    Constraints in the OpenMM simulation are not converted, and a warning will be issued
    if they are present.

    :param simulation: An OpenMM simulation object.

    :return: An ASE `MolecularDynamics` object representing the same system.
    """
    if simulation.system.getNumConstraints() > 0:
        logging.warning(CONSTRAINTS_UNSUPPORTED_MESSAGE)

    atoms = openmm_simulation_to_ase_atoms(simulation)

    timestep = simulation.integrator.getStepSize().value_in_unit(picoseconds)

    if isinstance(simulation.integrator, LangevinIntegrator):
        integrator: Callable[..., MolecularDynamics] = functools.partial(
            Langevin,
            timestep=timestep * _OpenMMToASE.time,
            temperature_K=simulation.integrator.getTemperature()._value,  # type: ignore
            friction=simulation.integrator.getFriction().value_in_unit(
                picoseconds ** (-1)
            )
            / _OpenMMToASE.time,
            fixcm=False,
        )
    elif isinstance(simulation.integrator, VerletIntegrator):
        logging.warning(INTEGRATOR_VERLET_MESSAGE)
        integrator = functools.partial(
            VelocityVerlet,
            timestep=timestep * _OpenMMToASE.time,
        )
    else:
        logging.warning(INTEGRATOR_NOT_LANGEVIN_MESSAGE)
        integrator = functools.partial(
            Langevin,
            timestep=timestep * _OpenMMToASE.time,
            temperature_K=300,
            friction=DEFAULT_LANGEVIN_FRICTION / _OpenMMToASE.time,
            fixcm=False,
        )

    return integrator(atoms=atoms)
