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
# Originally part of the narupa-ase package.
# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Modified under the terms of the GPL.

from threading import Lock
from typing import Any, List, Optional, Tuple

import numpy as np
from ase.atoms import Atoms
from ase.calculators.calculator import CalculatorSetupError, all_changes
from openmm.app import Simulation
from openmm.unit import angstrom

from narupatools.ase._units import UnitsASE
from narupatools.ase.calculators import Calculator
from narupatools.openmm import UnitsOpenMM

_OpenMMToASE = UnitsOpenMM >> UnitsASE


class OpenMMCalculator(Calculator):
    """
    Simple implementation of an ASE calculator for OpenMM.

    The context of the OpenMM simulation is used to compute forces and energies given a
    set of positions. When the ASE `Atoms` object has its positions changed by an
    integrator, these changes are pushed to the OpenMM context to enable the calculation
    of new forces and energies.
    """

    implemented_properties = ["energy", "forces"]

    ignored_changes = set(all_changes) - {"positions"}

    def __init__(self, simulation: Simulation, **kwargs: Any):
        """
        Create a calculator for the given simulation.

        :param simulation: OpenMM simulation to use as a calculator.
        :param atoms: Atoms object to which this calculator will be attached.
        :param kwargs: Dictionary of keywords to pass to the base ASE calculator.
        """
        super().__init__(**kwargs)
        self._simulation = simulation
        self._calculate_lock = Lock()

    def _calculate(  # noqa: D102
        self, atoms: Atoms, *, system_changes: List[str] = all_changes, **kwargs: Any
    ) -> None:
        if "positions" in system_changes:
            self._set_positions(atoms.positions)
            self.results["energy"], self.results["forces"] = self._calculate_openmm()

    def _calculate_openmm(self) -> Tuple[float, np.ndarray]:
        state = self._simulation.context.getState(getEnergy=True, getForces=True)
        energy = state.getPotentialEnergy()._value
        forces = state.getForces(asNumpy=True)._value
        return energy * _OpenMMToASE.energy, forces * _OpenMMToASE.force

    def _set_positions(self, positions: np.ndarray) -> None:
        self._simulation.context.setPositions(positions * angstrom)

    def assign_atoms(self, atoms: Optional[Atoms]) -> None:  # noqa: D102
        if (
            atoms is not None
            and len(atoms) != self._simulation.context.getSystem().getNumParticles()
        ):
            raise CalculatorSetupError(
                "OpenMM simulation and ASE Atoms have different sizes."
            )
