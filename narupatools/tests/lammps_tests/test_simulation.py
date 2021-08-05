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

import numpy as np
import pytest

lammps = pytest.importorskip("lammps")

from narupatools.lammps.simulation import (
    Compute,
    ComputeNotFoundError,
    InvalidComputeSpecificationError,
    LAMMPSSimulation,
    PropertyType,
    VariableStyle,
    VariableType,
)


@pytest.fixture(scope="module")
def simulation():
    return LAMMPSSimulation.from_file("./in.peptide", "./data.peptide")


@pytest.fixture(scope="module")
def simulation_dynamic():
    return LAMMPSSimulation.from_file("./in.peptide", "./data.peptide")


def test_extract_compute_missing(simulation):
    with pytest.raises(ComputeNotFoundError):
        simulation.extract_compute(
            "missing_key", VariableStyle.ATOM, VariableType.SCALAR
        )


def test_extract_compute_invalid(simulation):
    with pytest.raises(InvalidComputeSpecificationError):
        simulation.extract_compute("thermo_pe", VariableStyle.ATOM, VariableType.SCALAR)


def test_atom_property_valid(simulation):
    assert simulation.gather_atoms("f", PropertyType.DOUBLE, 3) is not None


def test_packages(simulation):
    assert isinstance(simulation.lammps_packages, list)
    assert "PYTHON" in simulation.lammps_packages


def test_computes(simulation):
    assert isinstance(simulation.computes, dict)
    assert (
        Compute(simulation=simulation, group="all", name="thermo_temp", style="temp")
        == simulation.computes["thermo_temp"]
    )
    assert (
        Compute(
            simulation=simulation, group="all", name="thermo_press", style="pressure"
        )
        == simulation.computes["thermo_press"]
    )
    assert (
        Compute(simulation=simulation, group="all", name="thermo_pe", style="pe")
        == simulation.computes["thermo_pe"]
    )


def test_energy(simulation):
    assert simulation.potential_energy == pytest.approx(-26662.02362458588)


def test_temperature(simulation):
    assert simulation.temperature == pytest.approx(275.0, abs=15.0)


def test_positions(simulation):
    assert len(simulation.positions) == 2004
    assert simulation.positions[0] == pytest.approx(
        np.array([4.399993, 5.852678, 3.67855])
    )


def test_velocities(simulation):
    assert len(simulation.velocities) == 2004


def test_forces(simulation):
    assert len(simulation.forces) == 2004


def test_masses(simulation):
    assert len(simulation.masses) == 2004
    assert simulation.masses[0] > 0


def test_set_positions(simulation_dynamic):
    positions = simulation_dynamic.positions.copy()
    positions[0] += [10.0, 0.0, 0.0]
    simulation_dynamic.positions = positions
    new_positions = simulation_dynamic.positions
    assert new_positions[0] == pytest.approx(positions[0])
