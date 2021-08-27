import contextlib
import re
import warnings
from typing import Generator

from lammps import OutputCapture

from narupatools.lammps import LAMMPSError
from narupatools.lammps.exceptions import AtomIDsNotDefinedError, UnknownPropertyNameError, LAMMPSWarning, \
    MissingInputScriptError, CannotOpenFileError, IllegalCommandError, UnknownCommandError, UnrecognizedStyleError

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


@contextlib.contextmanager
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