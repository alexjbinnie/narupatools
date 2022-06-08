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

"""IPython magick for using OpenMM."""
from typing import Any

from IPython.core.magic import Magics, line_magic, magics_class

from ._show import show


@magics_class
class NGLMagics(Magics):
    """IPython magic for use with NGL."""

    @line_magic
    def ngl(self, line: str) -> Any:
        """Display the given object using nglview."""
        widget = show(self.shell.ev(line))
        widget.add_unitcell()
        return widget
