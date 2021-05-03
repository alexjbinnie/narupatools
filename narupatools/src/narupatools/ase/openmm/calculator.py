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

"""ASE Calculator that interfaces with OpenMM to provide forces."""

from typing import Any, Collection, List, Optional, Tuple

import numpy as np
from ase.atoms import Atoms
from ase.calculators.calculator import Calculator, CalculatorSetupError, all_changes
from simtk.openmm.app import Simulation
from simtk.openmm.openmm import Context
from simtk.unit import angstrom

from narupatools.ase.calculators import CalculatorSetAtoms
from narupatools.ase.constraints.observer import ASEObserver
from narupatools.ase.units import UnitsASE
from narupatools.openmm.units import UnitsOpenMM
from narupatools.physics.typing import Vector3Array

_OpenMMToASE = UnitsOpenMM >> UnitsASE


class OpenMMCalculator(Calculator, CalculatorSetAtoms):
    """
    Simple implementation of a ASE calculator for OpenMM.

    The context of the OpenMM simulation is used to compute forces and energies given a
    set of positions. When the ASE `Atoms` object has its positions changed by an
    integrator, these changes are pushed to the OpenMM context to enable the calculation
    of new forces and energies.
    """

    implemented_properties = ["energy", "forces"]

    def __init__(self, simulation: Simulation, **kwargs: Any):
        """
        Create a calculator for the given simulation.

        :param simulation: OpenMM simulation to use as a calculator.
        :param kwargs: Dictionary of keywords to pass to the base ASE calculator.
        """
        super().__init__(**kwargs)
        self._context: Context = simulation.context
        self._positions_dirty: bool = True
        self._energy: Optional[float] = None
        self._forces: Optional[Vector3Array] = None
        self._position_observer = ASEObserver()
        self._position_observer.on_set_positions.add_callback(
            self._mark_positions_as_dirty
        )

    def set_atoms(self, /, atoms: Atoms) -> None:
        """
        Called when this is assigned using :ref:`Atoms.set_calculator`.

        :param atoms: ASE atoms object this calculator has been assigned to.
        """
        if self.atoms is not None:
            self.atoms.constraints.remove(self._position_observer)
        self.atoms = atoms
        self._ensure_atoms_valid(atoms)
        self.atoms.constraints.append(self._position_observer)
        self._mark_positions_as_dirty()

    def _mark_positions_as_dirty(self, **kwargs: Any) -> None:
        self._positions_dirty = True
        self.reset()

    def calculate(  # noqa: D102
        self,
        atoms: Optional[Atoms] = None,
        properties: Collection[str] = ("energy", "forces"),
        system_changes: List[str] = all_changes,
    ) -> None:
        if self.atoms is not None and (atoms is self.atoms or atoms is None):
            if self._positions_dirty:
                self._set_positions(self.atoms.positions)
                self._positions_dirty = False
                self._energy, self._forces = self._calculate_openmm()
        elif atoms is not None:
            self._ensure_atoms_valid(atoms)
            self._set_positions(atoms.positions)
            self._positions_dirty = True
            self._energy, self._forces = self._calculate_openmm()
            self.set_atoms(atoms)
        else:
            raise CalculatorSetupError("Atoms not set.")

        self.results["energy"] = self._energy
        self.results["forces"] = self._forces

    def _calculate_openmm(self) -> Tuple[float, np.ndarray]:
        state = self._context.getState(getEnergy=True, getForces=True)
        energy = state.getPotentialEnergy()._value
        forces = state.getForces(asNumpy=True)._value
        return energy * _OpenMMToASE.energy, forces * _OpenMMToASE.force

    def _set_positions(self, positions: np.ndarray) -> None:
        self._context.setPositions(positions * angstrom)

    def _ensure_atoms_valid(self, atoms: Atoms) -> None:
        if len(atoms) != self._context.getSystem().getNumParticles():
            raise CalculatorSetupError(
                "OpenMM simulation and ASE Atoms have different sizes."
            )
