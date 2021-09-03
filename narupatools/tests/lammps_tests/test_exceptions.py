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

from narupatools.lammps import LAMMPSSimulation
from narupatools.lammps._constants import VariableType
from narupatools.lammps.atom_properties import AtomProperty
from narupatools.lammps.exceptions import (
    CannotOpenFileError,
    IllegalCommandError,
    MissingInputScriptError,
    UnknownAtomPropertyError,
    UnknownCommandError,
    UnrecognizedStyleError,
)


def test_missing_input_file():
    simulation = LAMMPSSimulation.create_new("real")
    with pytest.raises(MissingInputScriptError):
        simulation.file("in.missing")


def test_present_input_file():
    simulation = LAMMPSSimulation.create_new("real")
    simulation.file("in.peptide")


def test_missing_data_file():
    simulation = LAMMPSSimulation.create_new("real")
    simulation.command("atom_style full")
    simulation.command("pair_style lj/charmm/coul/long 8.0 10.0 10.0")
    simulation.command("bond_style harmonic")
    simulation.command("angle_style charmm")
    simulation.command("dihedral_style charmm")
    simulation.command("improper_style harmonic")
    simulation.command("kspace_style pppm 0.0001")
    with pytest.raises(CannotOpenFileError):
        simulation.command("read_data data.missing")


def test_present_data_file():
    simulation = LAMMPSSimulation.create_new("real")
    simulation.command("atom_style full")
    simulation.command("pair_style lj/charmm/coul/long 8.0 10.0 10.0")
    simulation.command("bond_style harmonic")
    simulation.command("angle_style charmm")
    simulation.command("dihedral_style charmm")
    simulation.command("improper_style harmonic")
    simulation.command("kspace_style pppm 0.0001")
    simulation.command("read_data data.peptide")


def test_unrecognized_atom_style():
    simulation = LAMMPSSimulation.create_new("real")
    with pytest.raises(UnrecognizedStyleError):
        simulation.command("atom_style unknown")


def test_valid_atom_style():
    simulation = LAMMPSSimulation.create_new("real")
    simulation.command("atom_style full")


def test_unrecognized_pair_style():
    simulation = LAMMPSSimulation.create_new("real")
    with pytest.raises(UnrecognizedStyleError):
        simulation.command("pair_style unknown")


def test_valid_pair_style():
    simulation = LAMMPSSimulation.create_new("real")
    simulation.command("pair_style lj/charmm/coul/long 8.0 10.0 10.0")


def test_valid_bond_style():
    simulation = LAMMPSSimulation.create_new("real")
    simulation.command("atom_style full")
    simulation.command("pair_style lj/charmm/coul/long 8.0 10.0 10.0")
    simulation.command("bond_style harmonic")


def test_unknown_command():
    simulation = LAMMPSSimulation.create_new("real")
    with pytest.raises(UnknownCommandError):
        simulation.command("missing_command 2.0")


def test_illegal_command():
    simulation = LAMMPSSimulation.create_new("real")
    with pytest.raises(IllegalCommandError):
        simulation.command("atom_modify map blah")


def test_invalid_gather_atoms_name():
    simulation = LAMMPSSimulation.from_file("in.peptide")
    property = AtomProperty.define(
        "unknown_key", type=VariableType.DOUBLE, components=3
    )
    with pytest.raises(UnknownAtomPropertyError):
        simulation.gather_atoms(property)
