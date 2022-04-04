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

import re
import warnings

from narupatools.lammps import LAMMPSError
from narupatools.lammps.exceptions import (
    AtomIDsNotDefinedError,
    CannotOpenFileError,
    IllegalCommandError,
    LAMMPSWarning,
    MissingInputScriptError,
    UnknownCommandError,
    UnknownDataTypeError,
    UnrecognizedStyleError,
)

ILLEGAL_COMMAND_REGEX = re.compile(r"illegal \w+ command", re.IGNORECASE)
UNKNOWN_COMMAND_REGEX = re.compile(r"unknown command", re.IGNORECASE)
MISSING_INPUT_SCRIPT_REGEX = re.compile(
    r"cannot open input script [\w.]+: No such file or directory", re.IGNORECASE
)
CANNOT_OPEN_FILE_REGEX = re.compile(
    r"cannot open file [\w.]+: No such file or directory", re.IGNORECASE
)
UNRECOGNIZED_STYLE_REGEX = re.compile(r"Unrecognized \w+ style", re.IGNORECASE)


def _handle_error(message: str) -> None:
    if MISSING_INPUT_SCRIPT_REGEX.match(message) is not None:
        raise MissingInputScriptError(message)
    if CANNOT_OPEN_FILE_REGEX.match(message) is not None:
        raise CannotOpenFileError(message)
    if ILLEGAL_COMMAND_REGEX.match(message) is not None:
        raise IllegalCommandError(message)
    if UNKNOWN_COMMAND_REGEX.match(message) is not None:
        raise UnknownCommandError(message)
    if UNRECOGNIZED_STYLE_REGEX.match(message) is not None:
        raise UnrecognizedStyleError(message)
    raise LAMMPSError(message)


def handle_lammps_output_line(line: str) -> None:
    """Handle LAMMPS output and convert to errors if necessary."""
    if line.startswith("WARNING: "):
        warning = line[9:]
        if warning.startswith("Library error in lammps_gather_atoms"):
            raise AtomIDsNotDefinedError(func_name="gather_atoms")
        elif warning.startswith("Library error in lammps_scatter_atoms"):
            raise AtomIDsNotDefinedError(func_name="scatter_atoms")
        # elif warning.startswith("lammps_gather_atoms: unknown property name"):
        # raise UnknownPropertyNameError(warning)
        elif warning.startswith("lammps_gather_atoms: unsupported data type"):
            raise UnknownDataTypeError(warning)
        else:
            warnings.warn(LAMMPSWarning(warning))
    elif line.startswith("ERROR: "):
        error = line[7:]
        _handle_error(error)
