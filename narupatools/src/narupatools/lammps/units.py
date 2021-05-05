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

"""Unit conversions for LAMMPS."""

from typing import Dict

from narupatools.core.units import (
    UnitSystem,
    amu,
    angstrom,
    atmosphere,
    atomic_time_unit,
    atto,
    bar,
    bohr,
    calorie,
    centi,
    coulomb,
    debye,
    electronvolt,
    elementary_charge,
    femto,
    gram,
    hartree,
    kelvin,
    kilo,
    meter,
    micro,
    mole,
    nano,
    pascal,
    pico,
    poise,
    second,
    statcoulomb,
    volt,
)

UnitsLAMMPSReal = UnitSystem(
    length=angstrom,
    time=femto * second,
    mass=gram / mole,
    charge=elementary_charge,
    temperature=kelvin,
    pressure=atmosphere,
    energy=kilo * calorie / mole,
    force=kilo * calorie / (mole * angstrom),
    electric_field=volt / angstrom,
    density=gram / ((centi * meter) ** 3),
    density2d=gram / ((centi * meter) ** 2),
    dynamic_viscosity=poise,
)

UnitsLAMMPSMetal = UnitSystem(
    mass=gram / mole,
    length=angstrom,
    time=pico * second,
    energy=electronvolt,
    force=electronvolt / angstrom,
    temperature=kelvin,
    dynamic_viscosity=poise,
    electric_field=volt / angstrom,
    pressure=bar,
    charge=elementary_charge,
    density=gram / ((centi * meter) ** 3),
    density2d=gram / ((centi * meter) ** 2),
)

UnitsLAMMPSSI = UnitSystem(
    mass=kilo * gram, length=meter, time=second, temperature=kelvin, charge=coulomb
)

UnitsLAMMPSCGS = UnitSystem(
    mass=gram, length=centi * meter, time=second, temperature=kelvin, charge=statcoulomb
)

UnitsLAMMPSElectron = UnitSystem(
    mass=amu,
    length=bohr,
    time=femto * second,
    energy=hartree,
    charge=elementary_charge,
    temperature=kelvin,
    velocity=bohr / atomic_time_unit,
    force=hartree / bohr,
    pressure=pascal,
    dipole_moment=debye,
    electric_field=volt / (centi * meter),
)

UnitsLAMMPSMicro = UnitSystem(
    mass=pico * gram,
    length=micro * meter,
    time=micro * second,
    temperature=kelvin,
    charge=pico * coulomb,
    electric_field=volt / (micro * meter),
)

UnitsLAMMPSNano = UnitSystem(
    mass=atto * gram,
    length=nano * meter,
    time=nano * second,
    temperature=kelvin,
    charge=elementary_charge,
    electric_field=volt / (nano * meter),
)

_UNIT_SYSTEMS: Dict[str, UnitSystem] = {
    "real": UnitsLAMMPSReal,
    "metal": UnitsLAMMPSMetal,
    "si": UnitsLAMMPSSI,
    "cgs": UnitsLAMMPSCGS,
    "electron": UnitsLAMMPSElectron,
    "micro": UnitsLAMMPSMicro,
    "nano": UnitsLAMMPSNano,
}


def get_unit_system(system: str) -> UnitSystem:
    """
    Get the LAMMPS unit system for the given unit style.

    :param system: Name of the unit system as defined by LAMMPS.
    """
    return _UNIT_SYSTEMS[system]
