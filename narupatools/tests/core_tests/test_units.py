import math

import numpy as np
import pytest

from narupatools.ase._units import UnitsASE
from narupatools.core.units import (
    UnfinishedUnitOrQuantity,
    Unit,
    UnitsNarupa,
    amp,
    amu,
    angstrom,
    centi,
    coulomb,
    debye,
    degree,
    electronvolt,
    elementary_charge,
    erg,
    esu,
    gram,
    hartree,
    joule,
    kelvin,
    kilo,
    meter,
    milli,
    mole,
    nano,
    newton,
    pascal,
    pico,
    poise,
    radian,
    second,
    statcoulomb,
    volt,
)
from narupatools.mdanalysis._units import UnitsMDAnalysis
from narupatools.mdtraj._units import UnitsMDTraj
from narupatools.openmm._units import UnitsOpenMM


def test_general():
    assert 1.0 * newton == pytest.approx(1.0 * kilo * gram * meter / (second ** 2))
    assert 1.0 * pascal == pytest.approx(1.0 * newton / (meter ** 2))
    assert 1.0 * joule == pytest.approx(1.0 * newton * meter)

    assert 1.0 * hartree == pytest.approx(4.3597447222071e-18 * joule)

    assert 1.0 * volt == pytest.approx(1.0 * joule / coulomb)
    assert 1.0 * volt == pytest.approx(
        1.0 * kilo * gram * meter * meter / (second ** 3) / amp
    )

    assert 1.0 * erg == pytest.approx(
        1.0 * gram * ((centi * meter) ** 2) / (second ** 2)
    )
    assert 1.0 * erg == pytest.approx(1.0 * joule * 1e-7)

    assert 1.0 * statcoulomb == pytest.approx(1.0 * esu)
    assert 1.0 * statcoulomb == pytest.approx(3.33564095198e-10 * coulomb)

    assert 1.0 * poise == pytest.approx(
        0.1 * (meter ** -1) * kilo * gram * (second ** -1)
    )
    assert 1.0 * poise == pytest.approx(0.1 * pascal * second)

    assert 1.0 * debye == pytest.approx(1e-18 * statcoulomb * centi * meter)
    assert 1.0 * debye == pytest.approx(1e-10 * esu * angstrom)
    assert 1.0 * debye == pytest.approx(0.2081942 * elementary_charge * angstrom)
    assert 1.0 * debye == pytest.approx(3.33564095198e-30 * coulomb * meter)


def test_conversions():
    assert (meter >> angstrom) == pytest.approx(1e10)


def test_narupa():
    units = UnitsNarupa

    assert 1.0 * units.mass == pytest.approx(1.0)
    assert 1.0 * units.length == pytest.approx(1.0)
    assert 1.0 * units.time == pytest.approx(1.0)
    assert 1.0 * units.energy == pytest.approx(1.0)
    assert 1.0 * units.velocity == pytest.approx(1.0)
    assert 1.0 * units.force == pytest.approx(1.0)
    assert 1.0 * units.torque == pytest.approx(1.0)
    assert 1.0 * units.temperature == pytest.approx(1.0)
    assert 1.0 * units.pressure == pytest.approx(1.0)
    assert 1.0 * units.dynamic_viscosity == pytest.approx(1.0)
    assert 1.0 * units.charge == pytest.approx(1.0)
    assert 1.0 * units.dipole_moment == pytest.approx(1.0)
    assert 1.0 * units.electric_field == pytest.approx(1.0)
    assert 1.0 * units.density == pytest.approx(1.0)
    assert 1.0 * units.density2d == pytest.approx(1.0)

    assert 1.0 * kilo * joule / mole == pytest.approx(1.0)
    assert 1.0 * nano * meter == pytest.approx(1.0)
    assert 1.0 * pico * second == pytest.approx(1.0)
    assert 1.0 * elementary_charge == pytest.approx(1.0)
    assert 1.0 * kelvin == pytest.approx(1.0)

    # Check consistency of units
    assert 1.0 * units.energy == pytest.approx(
        1.0 * units.mass * (units.length ** 2) / (units.time ** 2)
    )


def test_ase():
    units = UnitsASE

    ase_time = angstrom * (amu / electronvolt) ** 0.5

    assert 1.0 * units.mass == pytest.approx(1.0 * amu)
    assert 1.0 * units.length == pytest.approx(1.0 * angstrom)
    assert 1.0 * units.time == pytest.approx(1.0 * ase_time)
    assert 1.0 * units.energy == pytest.approx(1.0 * electronvolt)
    assert 1.0 * units.velocity == pytest.approx(1.0 * angstrom / ase_time)
    assert 1.0 * units.force == pytest.approx(1.0 * amu * angstrom / (ase_time ** 2))
    assert 1.0 * units.force == pytest.approx(1.0 * electronvolt / angstrom)
    assert 1.0 * units.torque == pytest.approx(1.0 * electronvolt)
    assert 1.0 * units.temperature == pytest.approx(1.0 * kelvin)
    assert 1.0 * units.pressure == pytest.approx(1.0 * electronvolt / (angstrom ** 3))
    assert 1.0 * units.charge == pytest.approx(1.0 * elementary_charge)
    assert 1.0 * units.dipole_moment == pytest.approx(
        1.0 * elementary_charge * angstrom
    )
    assert 1.0 * units.density == pytest.approx(1.0 * amu / (angstrom ** 3))
    assert 1.0 * units.density2d == pytest.approx(1.0 * amu / (angstrom ** 2))

    # Check consistency of units
    assert 1.0 * units.energy == pytest.approx(
        1.0 * units.mass * (units.length ** 2) / (units.time ** 2)
    )


def test_mdtraj():
    units = UnitsMDTraj

    assert 1.0 * units.length == pytest.approx(1.0 * nano * meter)
    assert 1.0 * units.time == pytest.approx(1.0 * pico * second)
    assert 1.0 * units.angle == pytest.approx(1.0 * degree)
    assert 1.0 * units.energy == pytest.approx(1.0 * kilo * joule / mole)

    # Check consistency of units
    assert 1.0 * units.energy == pytest.approx(
        1.0 * units.mass * (units.length ** 2) / (units.time ** 2)
    )


def test_mdanalysis():
    units = UnitsMDAnalysis

    assert 1.0 * units.mass == pytest.approx(1.0 * amu)
    assert 1.0 * units.length == pytest.approx(1.0 * angstrom)
    assert 1.0 * units.time == pytest.approx(1.0 * pico * second)
    assert 1.0 * units.energy == pytest.approx(1.0 * kilo * joule / mole)
    assert 1.0 * units.force == pytest.approx(1.0 * kilo * joule / (mole * angstrom))
    assert 1.0 * units.angle == pytest.approx(1.0 * degree)
    assert 1.0 * units.velocity == pytest.approx(1.0 * angstrom / (pico * second))
    assert 1.0 * units.charge == pytest.approx(1.0 * elementary_charge)
    assert 1.0 * units.density == pytest.approx(1.0 * amu / (angstrom ** 3))

    # Check inconsistency of units
    assert 1.0 * units.energy != pytest.approx(
        1.0 * units.mass * (units.length ** 2) / (units.time ** 2)
    )


def test_openmm():
    units = UnitsOpenMM

    assert 1.0 * units.mass == pytest.approx(1.0 * nano * meter)
    assert 1.0 * units.time == pytest.approx(1.0 * pico * second)
    assert 1.0 * units.mass == pytest.approx(1.0 * amu)
    assert 1.0 * units.charge == pytest.approx(1.0 * elementary_charge)
    assert 1.0 * units.temperature == pytest.approx(1.0 * kelvin)
    assert 1.0 * units.angle == pytest.approx(1.0 * radian)
    assert 1.0 * units.energy == pytest.approx(1.0 * kilo * joule / mole)

    # Check consistency of units
    assert 1.0 * units.energy == pytest.approx(
        1.0 * units.mass * (units.length ** 2) / (units.time ** 2)
    )


@pytest.mark.parametrize("units", [UnitsOpenMM, UnitsMDAnalysis, UnitsMDTraj, UnitsASE])
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


@pytest.mark.parametrize("units", [UnitsOpenMM, UnitsMDAnalysis, UnitsMDTraj, UnitsASE])
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


@pytest.mark.parametrize("units", [UnitsOpenMM, UnitsMDAnalysis, UnitsMDTraj, UnitsASE])
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


def test_unit_divide():
    assert isinstance(meter / second, Unit)


def test_unit_multiply():
    assert isinstance(meter * second, Unit)


def test_unit_sqrt():
    with pytest.raises(TypeError):
        assert isinstance(math.sqrt(meter), Unit)


def test_unit_power():
    assert isinstance(meter ** 3, Unit)


def test_unit_inverse_power():
    assert isinstance(meter ** -2, Unit)


def test_unit_add():
    with pytest.raises(TypeError):
        assert isinstance(meter + second, Unit)


def test_unit_subtract():
    with pytest.raises(TypeError):
        assert isinstance(meter - second, Unit)


def test_unit_post_multiply_float():
    with pytest.raises(TypeError):
        _ = meter * 2.0


def test_unit_post_divide_float():
    with pytest.raises(TypeError):
        _ = meter / 2.0


def test_chain_prefixes():
    a = centi * milli * meter
    assert isinstance(a, Unit)
    b = 1.0 * a
    assert isinstance(b, float)


def test_unfinished():
    assert isinstance(2.0 * centi, UnfinishedUnitOrQuantity)
    assert isinstance(centi * meter * pico, UnfinishedUnitOrQuantity)
    assert isinstance(4.0 * milli * meter * centi, UnfinishedUnitOrQuantity)


def test_numpy_dtype():
    array = np.array([1.0, 2.0], dtype=float)
    assert (array * meter).dtype == float
    assert (array / meter).dtype == float
    assert (array * milli * meter).dtype == float
    assert (array / (centi * meter)).dtype == float
