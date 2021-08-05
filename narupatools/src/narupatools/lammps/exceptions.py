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

"""Custom exceptions that wrap LAMMPS output."""
from .constants import VariableStyle, VariableType


class LAMMPSError(RuntimeError):
    """Error raised by LAMMPS."""

    def __init__(self, message: str):
        super().__init__(message)


class CannotOpenFileError(LAMMPSError, FileNotFoundError):
    """Error raised when a file is not found by LAMMPS."""

    pass


class IllegalCommandError(LAMMPSError):
    """Error raised when a command is used incorrectly in LAMMPS."""

    pass


class UnknownPropertyNameError(LAMMPSError):
    """Error raised when lammps_gather_atoms encounters an unknown property name."""

    def __init__(self) -> None:
        super().__init__("Unknown property name when calling gather_atoms.")


class AtomIDsNotDefinedError(LAMMPSError):
    """Error raised when gather/scatter commands are used without atom ids."""

    def __init__(self, func_name: str):
        super().__init__(
            f"Atom IDs not defined or correctly ordered, so {func_name} cannot be "
            f"called. Use the 'atom_modify' command to allow this command to be used."
        )


class ComputeNotFoundError(LAMMPSError):
    """Error raised when a compute is requested which is not defined."""

    def __init__(self, key: str):
        super().__init__(f"No compute defined with key {key}")


class InvalidComputeSpecificationError(LAMMPSError):
    """Error raised when a compute is requested where the types are not defined."""

    def __init__(self, key: str, style: VariableStyle, type: VariableType):
        super().__init__(
            f"Invalid compute specification for {key}: Style {style} and type {type}"
        )


class UnknownAtomPropertyError(LAMMPSError):
    """Error raised when gather_atoms is used to get a property which does not exist."""

    def __init__(self, key: str):
        super().__init__(f"Attempted to get unknown atom property {key}.")
