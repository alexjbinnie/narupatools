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

from .units import joule, kelvin, coulomb, kilo, gram, meter, second

boltzmann_constant = 1.380649e-23 * joule / kelvin
"""Boltzmann constant in kilojoules per mole per kelvin."""

elementary_charge = 1.602176634e-19 * coulomb
"""Elementary charge constant, equal to 1 in Narupa units."""

electron_mass = 9.1093837015e-31 * kilo * gram
"""Mass of the electron in daltons."""

speed_of_light = 299792458 * meter / second
"""Speed of light in a vacuum in nanometers per picoseconds."""
