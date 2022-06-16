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

"""Technique for piping output, for use with LAMMPS."""

import os
import tempfile
from threading import Thread
from typing import Any


class PipedOutput:
    """Creates a temporary file that allows file writes to be piped into Python."""

    def __init__(self, handle_line: Any):
        self._filename = os.path.join(tempfile.mkdtemp(), "fifo")
        os.mkfifo(self.filename)  # type: ignore
        self._thread = Thread(target=self._thread_run, daemon=True)
        self._thread.start()
        self._handle_line = handle_line

    @property
    def filename(self) -> str:
        """Filename on the system that should be written to."""
        return self._filename

    def _thread_run(self) -> None:
        with open(self.filename, "r") as file:
            line = ""
            while True:
                data = file.read(1)
                if data == "\n":
                    self._handle_line(line)
                    line = ""
                else:
                    line += data

    def close(self) -> None:
        """Close the pipe."""
        os.unlink(self.filename)
