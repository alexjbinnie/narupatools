# ----------------------------------------------------------------------
#   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
#   https://www.lammps.org/ Sandia National Laboratories
#   Steve Plimpton, sjplimp@sandia.gov
#
#   Copyright (2003) Sandia Corporation.  Under the terms of Contract
#   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
#   certain rights in this software.  This software is distributed under
#   the GNU General Public License.
#
# -------------------------------------------------------------------------

################################################################################
# Alternative Python Wrapper
# Written by Richard Berger <richard.berger@temple.edu>
################################################################################

"""Captures output from stdout, for reading results from LAMMPS."""

from __future__ import annotations

import os
import select
from typing import Any


class OutputCapture:
    """Utility class to capture LAMMPS library output."""

    def __init__(self) -> None:
        self.stdout_pipe_read, self.stdout_pipe_write = os.pipe()
        self.stdout_fd = 1

    def __enter__(self) -> OutputCapture:
        self.stdout = os.dup(self.stdout_fd)
        os.dup2(self.stdout_pipe_write, self.stdout_fd)
        return self

    def __exit__(self, exc_type: Any, exc_value: Any, traceback: Any) -> None:
        os.dup2(self.stdout, self.stdout_fd)
        os.close(self.stdout)
        os.close(self.stdout_pipe_read)
        os.close(self.stdout_pipe_write)

    # check if we have more to read from the pipe
    def _more_data(self, pipe: int) -> bool:
        r, _, _ = select.select([pipe], [], [], 0)
        return bool(r)

    # read the whole pipe
    def _read_pipe(self, pipe: int) -> str:
        out = ""
        while self._more_data(pipe):
            out += os.read(pipe, 1024).decode()
        return out

    @property
    def output(self) -> str:
        """Output that has been written to stdout so far."""
        return self._read_pipe(self.stdout_pipe_read)
