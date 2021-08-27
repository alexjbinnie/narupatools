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

lammps = pytest.importorskip("lammps")

from narupatools.core.units import (
    Unit,
    UnitsNarupa,
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
    dyne,
    electronvolt,
    elementary_charge,
    erg,
    esu,
    femto,
    gram,
    hartree,
    joule,
    kelvin,
    kilo,
    meter,
    micro,
    mole,
    nano,
    newton,
    pascal,
    pico,
    poise,
    second,
    statcoulomb,
    statvolt,
    volt,
)
from narupatools.lammps._units import (
    UnitsLAMMPSCGS,
    UnitsLAMMPSElectron,
    UnitsLAMMPSMetal,
    UnitsLAMMPSMicro,
    UnitsLAMMPSNano,
    UnitsLAMMPSReal,
    UnitsLAMMPSSI,
)


def test_lammps_real():
    units = UnitsLAMMPSReal

    assert 1.0 * units.mass == pytest.approx(1.0 * gram / mole)
    assert 1.0 * units.length == pytest.approx(1.0 * angstrom)
    assert 1.0 * units.time == pytest.approx(1.0 * femto * second)
    assert 1.0 * units.energy == pytest.approx(1.0 * kilo * calorie / mole)
    assert 1.0 * units.velocity == pytest.approx(1.0 * angstrom / (femto * second))
    assert 1.0 * units.force == pytest.approx(1.0 * kilo * calorie / (mole * angstrom))
    assert 1.0 * units.torque == pytest.approx(1.0 * kilo * calorie / mole)
    assert 1.0 * units.temperature == pytest.approx(1.0 * kelvin)
    assert 1.0 * units.pressure == pytest.approx(1.0 * atmosphere)
    assert 1.0 * units.dynamic_viscosity == pytest.approx(1.0 * poise)
    assert 1.0 * units.charge == pytest.approx(1.0 * elementary_charge)
    assert 1.0 * units.dipole_moment == pytest.approx(
        1.0 * elementary_charge * angstrom
    )
    assert 1.0 * units.electric_field == pytest.approx(1.0 * volt / angstrom)
    assert 1.0 * units.density == pytest.approx(1.0 * gram / ((centi * meter) ** 3))
    assert 1.0 * units.density2d == pytest.approx(1.0 * gram / ((centi * meter) ** 2))

    # Check inconsistency of units
    assert 1.0 * units.energy != pytest.approx(
        1.0 * units.mass * (units.length ** 2) / (units.time ** 2)
    )


def test_lammps_metal():
    units = UnitsLAMMPSMetal

    assert 1.0 * units.mass == pytest.approx(1.0 * gram / mole)
    assert 1.0 * units.length == pytest.approx(1.0 * angstrom)
    assert 1.0 * units.time == pytest.approx(1.0 * pico * second)
    assert 1.0 * units.energy == pytest.approx(1.0 * electronvolt)
    assert 1.0 * units.velocity == pytest.approx(1.0 * angstrom / (pico * second))
    assert 1.0 * units.force == pytest.approx(1.0 * electronvolt / angstrom)
    assert 1.0 * units.torque == pytest.approx(1.0 * electronvolt)
    assert 1.0 * units.temperature == pytest.approx(1.0 * kelvin)
    assert 1.0 * units.pressure == pytest.approx(1.0 * bar)
    assert 1.0 * units.dynamic_viscosity == pytest.approx(1.0 * poise)
    assert 1.0 * units.charge == pytest.approx(1.0 * elementary_charge)
    assert 1.0 * units.dipole_moment == pytest.approx(
        1.0 * elementary_charge * angstrom
    )
    assert 1.0 * units.electric_field == pytest.approx(1.0 * volt / angstrom)
    assert 1.0 * units.density == pytest.approx(1.0 * gram / ((centi * meter) ** 3))
    assert 1.0 * units.density2d == pytest.approx(1.0 * gram / ((centi * meter) ** 2))

    # Check inconsistency of units
    assert 1.0 * units.energy != pytest.approx(
        1.0 * units.mass * (units.length ** 2) / (units.time ** 2)
    )


def test_lammps_si():
    units = UnitsLAMMPSSI

    assert 1.0 * units.mass == pytest.approx(1.0 * kilo * gram)
    assert 1.0 * units.length == pytest.approx(1.0 * meter)
    assert 1.0 * units.time == pytest.approx(1.0 * second)
    assert 1.0 * units.energy == pytest.approx(1.0 * joule)
    assert 1.0 * units.velocity == pytest.approx(1.0 * meter / second)
    assert 1.0 * units.force == pytest.approx(1.0 * newton)
    assert 1.0 * units.torque == pytest.approx(1.0 * newton * meter)
    assert 1.0 * units.temperature == pytest.approx(1.0 * kelvin)
    assert 1.0 * units.pressure == pytest.approx(1.0 * pascal)
    assert 1.0 * units.dynamic_viscosity == pytest.approx(1.0 * pascal * second)
    assert 1.0 * units.charge == pytest.approx(1.0 * coulomb)
    assert 1.0 * units.dipole_moment == pytest.approx(1.0 * coulomb * meter)
    assert 1.0 * units.electric_field == pytest.approx(1.0 * volt / meter)
    assert 1.0 * units.density == pytest.approx(1.0 * kilo * gram / (meter ** 3))
    assert 1.0 * units.density2d == pytest.approx(1.0 * kilo * gram / (meter ** 2))

    # Check consistency of units
    assert 1.0 * units.energy == pytest.approx(
        1.0 * units.mass * (units.length ** 2) / (units.time ** 2)
    )


def test_lammps_cgs():
    units = UnitsLAMMPSCGS

    assert 1.0 * units.mass == pytest.approx(1.0 * gram)
    assert 1.0 * units.length == pytest.approx(1.0 * centi * meter)
    assert 1.0 * units.time == pytest.approx(1.0 * second)
    assert 1.0 * units.energy == pytest.approx(1.0 * erg)
    assert 1.0 * units.velocity == pytest.approx(1.0 * centi * meter / second)
    assert 1.0 * units.force == pytest.approx(1.0 * dyne)
    assert 1.0 * units.torque == pytest.approx(1.0 * dyne * centi * meter)
    assert 1.0 * units.temperature == pytest.approx(1.0 * kelvin)
    assert 1.0 * units.pressure == pytest.approx(1.0 * dyne / ((centi * meter) ** 2))
    assert 1.0 * units.pressure == pytest.approx(1e-6 * bar)
    assert 1.0 * units.dynamic_viscosity == pytest.approx(1.0 * poise)
    assert 1.0 * units.charge == pytest.approx(1.0 * statcoulomb)
    assert 1.0 * units.dipole_moment == pytest.approx(1.0 * statcoulomb * centi * meter)
    assert 1.0 * units.dipole_moment == pytest.approx(1e18 * debye)
    assert 1.0 * units.electric_field == pytest.approx(1.0 * statvolt / (centi * meter))
    assert 1.0 * units.electric_field == pytest.approx(1.0 * dyne / esu)
    assert 1.0 * units.density == pytest.approx(1.0 * gram / ((centi * meter) ** 3))
    assert 1.0 * units.density2d == pytest.approx(1.0 * gram / ((centi * meter) ** 2))

    # Check consistency of units
    assert 1.0 * units.energy == pytest.approx(
        1.0 * units.mass * (units.length ** 2) / (units.time ** 2)
    )


def test_lammps_electron():
    units = UnitsLAMMPSElectron

    assert 1.0 * units.mass == pytest.approx(1.0 * amu)
    assert 1.0 * units.length == pytest.approx(1.0 * bohr)
    assert 1.0 * units.time == pytest.approx(1.0 * femto * second)
    assert 1.0 * units.energy == pytest.approx(1.0 * hartree)
    assert 1.0 * units.velocity == pytest.approx(1.0 * bohr / atomic_time_unit)
    assert 1.0 * units.force == pytest.approx(1.0 * hartree / bohr)
    assert 1.0 * units.temperature == pytest.approx(1.0 * kelvin)
    assert 1.0 * units.pressure == pytest.approx(1.0 * pascal)
    assert 1.0 * units.charge == pytest.approx(1.0 * elementary_charge)
    assert 1.0 * units.dipole_moment == pytest.approx(1.0 * debye)
    assert 1.0 * units.electric_field == pytest.approx(1.0 * volt / (centi * meter))

    # Check inconsistency of units
    assert 1.0 * units.energy != pytest.approx(
        1.0 * units.mass * (units.length ** 2) / (units.time ** 2)
    )


def test_lammps_micro():
    units = UnitsLAMMPSMicro

    assert 1.0 * units.mass == pytest.approx(1.0 * pico * gram)
    assert 1.0 * units.length == pytest.approx(1.0 * micro * meter)
    assert 1.0 * units.time == pytest.approx(1.0 * micro * second)
    assert 1.0 * units.energy == pytest.approx(
        1.0 * pico * gram * ((micro * meter) ** 2) / ((micro * second) ** 2)
    )
    assert 1.0 * units.velocity == pytest.approx(
        1.0 * (micro * meter) / (micro * second)
    )
    assert 1.0 * units.force == pytest.approx(
        1.0 * pico * gram * (micro * meter) / ((micro * second) ** 2)
    )
    assert 1.0 * units.torque == pytest.approx(
        1.0 * pico * gram * ((micro * meter) ** 2) / ((micro * second) ** 2)
    )
    assert 1.0 * units.temperature == pytest.approx(1.0 * kelvin)
    assert 1.0 * units.pressure == pytest.approx(
        1.0 * pico * gram / (micro * meter * ((micro * second) ** 2))
    )
    assert 1.0 * units.dynamic_viscosity == pytest.approx(
        1.0 * (pico * gram) / (micro * meter * micro * second)
    )
    assert 1.0 * units.charge == pytest.approx(1.0 * pico * coulomb)
    assert 1.0 * units.dipole_moment == pytest.approx(
        1.0 * pico * coulomb * micro * meter
    )
    assert 1.0 * units.electric_field == pytest.approx(1.0 * volt / (micro * meter))
    assert 1.0 * units.density == pytest.approx(
        1.0 * pico * gram / ((micro * meter) ** 3)
    )
    assert 1.0 * units.density2d == pytest.approx(
        1.0 * pico * gram / ((micro * meter) ** 2)
    )

    # Check consistency of units
    assert 1.0 * units.energy == pytest.approx(
        1.0 * units.mass * (units.length ** 2) / (units.time ** 2)
    )


def test_lammps_nano():
    units = UnitsLAMMPSNano

    assert 1.0 * units.mass == pytest.approx(1.0 * atto * gram)
    assert 1.0 * units.length == pytest.approx(1.0 * nano * meter)
    assert 1.0 * units.time == pytest.approx(1.0 * nano * second)
    assert 1.0 * units.energy == pytest.approx(
        1.0 * atto * gram * ((nano * meter) ** 2) / ((nano * second) ** 2)
    )
    assert 1.0 * units.velocity == pytest.approx(1.0 * (nano * meter) / (nano * second))
    assert 1.0 * units.force == pytest.approx(
        1.0 * atto * gram * (nano * meter) / ((nano * second) ** 2)
    )
    assert 1.0 * units.torque == pytest.approx(
        1.0 * atto * gram * ((nano * meter) ** 2) / ((nano * second) ** 2)
    )
    assert 1.0 * units.temperature == pytest.approx(1.0 * kelvin)
    assert 1.0 * units.pressure == pytest.approx(
        1.0 * atto * gram / (nano * meter * ((nano * second) ** 2))
    )
    assert 1.0 * units.dynamic_viscosity == pytest.approx(
        1.0 * (atto * gram) / (nano * meter * nano * second)
    )
    assert 1.0 * units.charge == pytest.approx(1.0 * elementary_charge)
    assert 1.0 * units.dipole_moment == pytest.approx(
        1.0 * elementary_charge * nano * meter
    )
    assert 1.0 * units.electric_field == pytest.approx(1.0 * volt / (nano * meter))
    assert 1.0 * units.density == pytest.approx(
        1.0 * atto * gram / ((nano * meter) ** 3)
    )
    assert 1.0 * units.density2d == pytest.approx(
        1.0 * atto * gram / ((nano * meter) ** 2)
    )

    # Check consistency of units
    assert 1.0 * units.energy == pytest.approx(
        1.0 * units.mass * (units.length ** 2) / (units.time ** 2)
    )


@pytest.mark.parametrize(
    "units",
    [
        UnitsLAMMPSSI,
        UnitsLAMMPSNano,
        UnitsLAMMPSMetal,
        UnitsLAMMPSReal,
        UnitsLAMMPSCGS,
        UnitsLAMMPSElectron,
        UnitsLAMMPSMicro,
    ],
)
@pytest.mark.parametrize(
    "quantity",
    [
        "mass",
        "length",
        "time",
        "energy",
        "velocity",
        "force",
        "torque",
        "temperature",
        "pressure",
        "dynamic_viscosity",
        "charge",
        "dipole_moment",
        "density",
        "density2d",
    ],
)
def test_conversion_to_narupa(units, quantity):
    conversion = units >> UnitsNarupa
    conversion2 = UnitsNarupa << units
    assert 1.0 * getattr(conversion, quantity) == pytest.approx(
        1.0 * getattr(units, quantity)
    )
    assert 1.0 * getattr(conversion2, quantity) == pytest.approx(
        1.0 * getattr(units, quantity)
    )


@pytest.mark.parametrize(
    "units",
    [
        UnitsLAMMPSSI,
        UnitsLAMMPSNano,
        UnitsLAMMPSMetal,
        UnitsLAMMPSReal,
        UnitsLAMMPSCGS,
        UnitsLAMMPSElectron,
        UnitsLAMMPSMicro,
    ],
)
@pytest.mark.parametrize(
    "quantity",
    [
        "mass",
        "length",
        "time",
        "energy",
        "velocity",
        "force",
        "torque",
        "temperature",
        "pressure",
        "dynamic_viscosity",
        "charge",
        "dipole_moment",
        "density",
        "density2d",
    ],
)
def test_conversion_from_narupa(units, quantity):
    conversion = UnitsNarupa >> units
    conversion2 = units << UnitsNarupa
    assert 1.0 * getattr(conversion, quantity) == pytest.approx(
        1.0 / getattr(units, quantity)
    )
    assert 1.0 * getattr(conversion2, quantity) == pytest.approx(
        1.0 / getattr(units, quantity)
    )


@pytest.mark.parametrize(
    "units",
    [
        UnitsLAMMPSSI,
        UnitsLAMMPSNano,
        UnitsLAMMPSMetal,
        UnitsLAMMPSReal,
        UnitsLAMMPSCGS,
        UnitsLAMMPSElectron,
        UnitsLAMMPSMicro,
    ],
)
@pytest.mark.parametrize(
    "quantity",
    [
        "mass",
        "length",
        "time",
        "energy",
        "velocity",
        "force",
        "torque",
        "temperature",
        "pressure",
        "dynamic_viscosity",
        "charge",
        "dipole_moment",
        "density",
        "density2d",
    ],
)
def test_unit_system_is_unit(units, quantity):
    assert isinstance(getattr(units, quantity), Unit)
