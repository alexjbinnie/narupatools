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
import re
import warnings
from contextlib import contextmanager
from typing import Generator

from lammps import OutputCapture

from ._constants import VariableDimension, VariableStyle

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


@contextmanager
def catch_lammps_warnings_and_exceptions() -> Generator[None, None, None]:
    """
    Capture output and raises logged warnings and errors in a pythonic way.

    Any line starting with 'WARNING: ' will be raised as a LAMMPSWarning, except certain
    warnings which are instead treated as errors.

    Any line starting with 'ERROR:' will be raised as a LAMMPSError, unless a more
    specific error message exists.

    :raises AtomIDsNotDefinedError: Library error raises in gather/scatter atoms.
    :raises LAMMPSError: An ERROR: is logged to the console by LAMMPS.
    """
    with OutputCapture() as o:
        try:
            yield
        except Exception as e:
            if isinstance(e, LAMMPSError):
                raise e
            if e.args[0].startswith("ERROR on proc 0: "):
                _handle_error(e.args[0][17:])
            if e.args[0].startswith("ERROR: "):
                _handle_error(e.args[0][7:])
            raise LAMMPSError(e.args[0])
        output = o.output
    for line in output.splitlines():
        if line.startswith("WARNING: "):
            warning = line[9:]
            if warning.startswith("Library error in lammps_gather_atoms"):
                raise AtomIDsNotDefinedError(func_name="gather_atoms")
            elif warning.startswith("Library error in lammps_scatter_atoms"):
                raise AtomIDsNotDefinedError(func_name="scatter_atoms")
            elif warning.startswith("lammps_gather_atoms: unknown property name"):
                raise UnknownPropertyNameError(warning)
            else:
                warnings.warn(LAMMPSWarning(warning))
        elif line.startswith("ERROR: "):
            error = line[7:]
            _handle_error(error)


class UntestedVersionWarning(UserWarning):
    """Warning raised when a version of LAMMPS is used that may not be supported."""

    pass


class LAMMPSWarning(UserWarning):
    """Warning raised by LAMMPS."""

    pass


class LAMMPSError(RuntimeError):
    """Error raised by LAMMPS."""

    pass


class CannotOpenFileError(LAMMPSError, FileNotFoundError):
    """Error raised when a file is not found by LAMMPS."""

    pass


class MissingInputScriptError(CannotOpenFileError):
    """Error raised when a file is not found by LAMMPS."""

    pass


class UnknownCommandError(LAMMPSError):
    """Error raised when an unknown command is used in LAMMPS."""

    pass


class IllegalCommandError(LAMMPSError):
    """Error raised when a command is used incorrectly in LAMMPS."""

    pass


class UnrecognizedStyleError(LAMMPSError):
    """Error raised when a style not recognized by LAMMPS is used."""

    pass


class UnknownPropertyNameError(LAMMPSError):
    """Error raised when lammps_gather_atoms encounters an unknown property name."""


class AtomIDsNotDefinedError(LAMMPSError):
    """Error raised when gather/scatter commands are used without atom ids."""

    def __init__(self, func_name: str):
        super().__init__(
            f"Atom IDs not defined or correctly ordered, so {func_name} cannot be "
            f"called. Use the 'atom_modify' command to allow this command to be used."
        )


class VariableNotFoundError(LAMMPSError):
    """Error raised when a variable is requested which is not defined."""

    def __init__(self, key: str):
        super().__init__(f"No variable defined with key {key}")


class SettingNotFoundError(LAMMPSError):
    """Error raised when a setting is requested which is not defined."""

    def __init__(self, key: str):
        super().__init__(f"No setting defined with key {key}")


class GlobalNotFoundError(LAMMPSError):
    """Error raised when a global is requested which is not defined."""

    def __init__(self, key: str):
        super().__init__(f"No global defined with key {key}")


class InvalidThermoKeywordError(LAMMPSError):
    """Error raised when an invalid thermo keyword is specified."""

    pass


class ComputeNotFoundError(LAMMPSError):
    """Error raised when a compute is requested which is not defined."""

    def __init__(self, key: str):
        super().__init__(f"No compute defined with key {key}")


class FixNotFoundError(LAMMPSError):
    """Error raised when a fix is requested which is not defined."""

    def __init__(self, key: str):
        super().__init__(f"No fix defined with key {key}")


class InvalidComputeSpecificationError(LAMMPSError):
    """Error raised when a compute is requested where the types are not defined."""

    def __init__(self, key: str, style: VariableStyle, type: VariableDimension):
        super().__init__(
            f"Invalid compute specification for {key}: Style {style} and type {type}"
        )


class UnknownAtomPropertyError(LAMMPSError):
    """Error raised when gather_atoms is used to get a property which does not exist."""

    def __init__(self, key: str):
        super().__init__(f"Attempted to get unknown atom property {key}.")
