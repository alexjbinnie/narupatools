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

pytest.importorskip("lammps")


def test_initial_computes(lammps):
    assert len(lammps.computes) == 4
    assert set(lammps.computes) == {
        "thermo_temp",
        "thermo_press",
        "thermo_pe",
        "1_temp",
    }


def test_change_computes(lammps):
    lammps.command("compute 1 all ke")
    assert len(lammps.computes) == 5
    assert "1" in lammps.computes
    lammps.command("uncompute 1")
    assert len(lammps.computes) == 4
    assert "1" not in lammps.computes


def test_initial_dumps(lammps):
    assert len(lammps.dumps) == 0


def test_change_dumps(lammps):
    lammps.command("dump 1 all atom 100 dump.atom")
    assert len(lammps.dumps) == 1
    assert "1" in lammps.dumps
    lammps.command("undump 1")
    assert len(lammps.dumps) == 0
    assert "1" not in lammps.dumps


def test_initial_fixes(lammps):
    assert len(lammps.fixes) == 2
    assert set(lammps.fixes) == {"1", "2"}


def test_change_fixes(lammps):
    lammps.command("fix 3 all momentum 1 linear 1 1 0")
    assert len(lammps.fixes) == 3
    assert "3" in lammps.fixes
    lammps.command("unfix 3")
    assert len(lammps.fixes) == 2
    assert "3" not in lammps.fixes


def test_initial_groups(lammps):
    assert len(lammps.groups) == 1
    assert set(lammps.groups) == {"all"}


def test_change_groups(lammps):
    lammps.command("group water type 3 4")
    assert len(lammps.groups) == 2
    assert "water" in lammps.groups
    lammps.command("group water delete")
    assert len(lammps.groups) == 1
    assert "water" not in lammps.groups


def test_initial_molecules(lammps):
    assert len(lammps.molecules) == 0


def test_initial_regions(lammps):
    assert len(lammps.regions) == 0


def test_change_regions(lammps):
    lammps.command("region my_region block -3.0 5.0 INF 10.0 INF INF")
    assert len(lammps.regions) == 1
    assert "my_region" in lammps.regions
    lammps.command("region my_region delete")
    assert len(lammps.regions) == 0
    assert "my_region" not in lammps.regions


def test_initial_variables(lammps):
    assert len(lammps.variables) == 0


def test_change_variables(lammps):
    lammps.command("variable beta equal temp/3.0")
    assert len(lammps.variables) == 1
    assert "beta" in lammps.variables
    lammps.command("variable beta delete")
    assert len(lammps.variables) == 0
    assert "beta" not in lammps.variables
