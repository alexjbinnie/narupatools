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

"""ASE calculator that uses LAMMPS to calculate energies and forces."""

from threading import Lock
from typing import Any, Collection, List, Optional

import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes

from narupatools.ase import UnitsASE
from narupatools.ase.constraints import ASEObserver
from narupatools.physics.typing import Vector3Array
from narupatools.physics.units import UnitsNarupa

from ._simulation import LAMMPSSimulation

_NarupaToASE = UnitsNarupa >> UnitsASE
_ASEToNarupa = UnitsASE >> UnitsNarupa


class LAMMPSCalculator(Calculator):
    """
    ASE Calculator that interfaces with a given LAMMPS simulation.

    When forces or energies are requested, positions and velocities are sent to the
    LAMMPS simulation.
    """

    implemented_properties = ["energy", "forces"]

    def __init__(
        self, simulation: LAMMPSSimulation, atoms: Optional[Atoms] = None, **kwargs: Any
    ):
        super().__init__(**kwargs)
        self._simulation = simulation
        self._lammps_lock = Lock()
        self._atoms = atoms

        if atoms is not None:
            _position_observer = ASEObserver.get_or_create(atoms)
            _position_observer.on_set_positions.add_callback(
                self._mark_positions_as_dirty
            )
        self._positions_dirty = True

        self._energy = 0.0
        if self._atoms is None:
            self._forces = np.zeros((3, 0))
        else:
            self._forces = np.zeros((3, len(self._atoms)))

    @property
    def simulation(self) -> LAMMPSSimulation:
        """LAMMPS Simulation that provides energies and forces for this calculator."""
        return self._simulation

    def calculate(  # noqa: D102
        self,
        atoms: Optional[Atoms] = None,
        properties: Collection[str] = ("energy", "forces"),
        system_changes: List[str] = all_changes,
    ) -> None:
        if self._atoms is not None and (atoms is self._atoms or atoms is None):
            if self._positions_dirty:
                with self._lammps_lock:
                    self._set_positions(self._atoms.get_positions())
                    self._set_velocities(self._atoms.get_velocities())
                    self._run_system()
                    self._energy = self._extract_potential_energy()
                    self._forces = self._extract_forces()
                self._positions_dirty = False
        elif atoms is not None:
            with self._lammps_lock:
                self._set_positions(atoms.get_positions())
                self._run_system()
                self._energy = self._extract_potential_energy()
                self._forces = self._extract_forces()
            self._positions_dirty = True
        else:
            raise ValueError("No atoms provided")

        self.results["energy"] = self._energy
        self.results["forces"] = self._forces

    def _mark_positions_as_dirty(self, **kwargs: Any) -> None:
        self._positions_dirty = True

    def _extract_potential_energy(self) -> float:
        return self._simulation.potential_energy * _NarupaToASE.energy

    def _set_positions(self, positions: Vector3Array) -> None:
        self._simulation.positions = positions * _ASEToNarupa.length

    def _set_velocities(self, velocities: Vector3Array) -> None:
        self._simulation.velocities = velocities * _ASEToNarupa.velocity

    def _run_system(self) -> None:
        self._simulation.run(0)

    def _extract_forces(self) -> Vector3Array:
        return self._simulation.forces * _NarupaToASE.force
