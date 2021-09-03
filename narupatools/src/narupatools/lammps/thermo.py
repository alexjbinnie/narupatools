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

"""Thermo keywords defined by LAMMPS."""

from dataclasses import dataclass

from narupatools.lammps._wrapper import Extractable, LAMMPSWrapper, _TReturnType


@dataclass
class FloatThermoKeyword(Extractable[float]):
    """Thermo keyword that returns a float."""

    def extract(self, lammps: LAMMPSWrapper) -> float:
        return lammps.get_thermo(self.key)

    key: str


@dataclass
class IntegerThermoKeyword(Extractable[int]):
    """Thermo keyword that returns an integer."""

    def extract(self, lammps: LAMMPSWrapper) -> int:
        return int(lammps.get_thermo(self.key))

    key: str


Step = IntegerThermoKeyword("step")
"""Current timestep, or iteration count when a minimization is being performed."""

Elapsed = IntegerThermoKeyword("elapsed")
"""Number of timesteps elapsed since the beginning of this run."""

ElapsedLong = IntegerThermoKeyword("elaplong")
"""Number of timesteps elapsed since the beginning of an initial run in a series of runs."""

TimestepSize = FloatThermoKeyword("dt")
"""Current timestep size in LAMMPS internal time units."""

Time = FloatThermoKeyword("time")
"""Current elapsed simulation time in LAMMPS internal time units."""

CPUTime = FloatThermoKeyword("cpu")
"""Elapsed CPU seconds since the beginning of this run."""

SimulationTimePerCPUSecond = FloatThermoKeyword("tpcpu")
"""Simulation time (in LAMMPS internal time units) per CPU second."""

TimestepsPerCPUSecond = FloatThermoKeyword("spcpu")
"""Number of timesteps per CPU second."""

CPURemain = FloatThermoKeyword("cpuremain")
"""Estimate of the CPU time remaining in the current run, based on the time elapsed thus far."""

Partition = IntegerThermoKeyword("part")
"""Index of the partition this output and this file corresponds to,."""

TimeRemaining = FloatThermoKeyword("timeremain")
"""Seconds remaining when a timeout has been configured via the `timer timeout` command."""

NumberOfAtoms = IntegerThermoKeyword("atoms")
"""Number of atoms in the simulation."""

NumberOfBonds = IntegerThermoKeyword("bonds")
"""Number of bonds in the simulation."""

NumberOfAngles = IntegerThermoKeyword("angles")
"""Number of angles in the simulation."""

NumberOfDihedrals = IntegerThermoKeyword("dihedrals")
"""Number of dihedrals in the simulation."""

NumberOfImpropers = IntegerThermoKeyword("impropers")
"""Number of impropers in the simulation."""

Temperature = FloatThermoKeyword("temperature")
"""Current temperature of the simulation."""

Pressure = FloatThermoKeyword("press")
"""Current pressure of the simulation."""

PotentialEnergy = FloatThermoKeyword("pe")
"""Current potential energy of the simulation."""

KineticEnergy = FloatThermoKeyword("ke")
"""Current kinetic energy of the simulation."""

TotalEnergy = FloatThermoKeyword("etotal")
"""Total energy of the simulation."""

Enthalpy = FloatThermoKeyword("enthalpy")
"""Enthalpy of the simulation."""

VanDerWaalsEnergy = FloatThermoKeyword("evdwl")
"""Van der Waals pairwise energy."""

CoulombicEnergy = FloatThermoKeyword("ecoul")
"""Coulombic pairwise energy."""

PairwiseEnergy = FloatThermoKeyword("epair")
"""Pairwise energy."""

BondEnergy = FloatThermoKeyword("ebond")
"""Bond energy."""

AngleEnergy = FloatThermoKeyword("eangle")
"""Angle energy."""

DihedralEnergy = FloatThermoKeyword("edihed")
"""Dihedral energy."""

ImproperEnergy = FloatThermoKeyword("eimp")
"""Improper energy."""

MolecularEnergy = FloatThermoKeyword("emol")
"""Molecular energy."""

LongRangeKSpaceEnergy = FloatThermoKeyword("elong")
"""Long-range k-space energy."""

VanDerWaalsTailEnergy = FloatThermoKeyword("etail")
"""Van der Waals energy long-range tail correction."""

Volume = FloatThermoKeyword("vol")
"""Volume of the simulation."""

Density = FloatThermoKeyword("density")
"""Density of the simulation."""
