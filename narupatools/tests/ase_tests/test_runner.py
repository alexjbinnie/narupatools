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
from ase.calculators.calculator import Calculator
from ase.md import Langevin

from narupatools.ase import ASEDynamics, UnitsASE
from narupatools.ase import NullCalculator
from narupatools.core import UnitsNarupa

ASEToNarupa = UnitsASE >> UnitsNarupa


def test_from_ase_atoms_sets_null_calculator(ethane_atoms):
    ase = ASEDynamics.create_langevin(ethane_atoms)
    assert isinstance(ase.atoms.calc, NullCalculator)


class CustomCalculator(Calculator):
    pass


def test_from_ase_atoms_leaves_existing_calculator(ethane_atoms):
    calculator = CustomCalculator()
    ethane_atoms.set_calculator(calculator)
    ase = ASEDynamics.create_langevin(ethane_atoms)
    assert ase.atoms.calc == calculator


def test_from_ase_dynamics(villin_ase_langevin_dynamics, villin_mda_universe_readonly):
    dynamics = ASEDynamics(
        villin_ase_langevin_dynamics, universe=villin_mda_universe_readonly
    )
    assert dynamics.molecular_dynamics == villin_ase_langevin_dynamics


@pytest.mark.parametrize("friction", [1e-4, 1e-2])
@pytest.mark.parametrize("timestep", [1e-3, 1e-1])
@pytest.mark.parametrize("temperature", [300, 100])
def test_from_ase_langevin(villin_ase_atoms, friction, timestep, temperature):
    dynamics = ASEDynamics.create_langevin(
        villin_ase_atoms, friction=friction, temperature=temperature, timestep=timestep
    )
    assert isinstance(dynamics.molecular_dynamics, Langevin)
    assert dynamics.molecular_dynamics.todict()["temperature_K"] == pytest.approx(
        temperature
    )
    assert dynamics.molecular_dynamics.todict()[
        "friction"
    ] / ASEToNarupa.time == pytest.approx(friction)
    assert dynamics.molecular_dynamics.todict()[
        "timestep"
    ] * ASEToNarupa.time == pytest.approx(timestep)
