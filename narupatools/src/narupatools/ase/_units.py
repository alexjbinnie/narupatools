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

"""Unit conversions for ASE."""

from narupatools.core.units import (
    UnitSystem,
    amu,
    angstrom,
    electronvolt,
    elementary_charge,
    kelvin,
)

UnitsASE = UnitSystem(
    length=angstrom,
    mass=amu,
    temperature=kelvin,
    time=angstrom * (amu / electronvolt) ** 0.5,
    charge=elementary_charge,
)
r"""
Unit system for the ASE library.

ASE uses a set of consistent units. The defining units of ASE are the angstrom (length),
amu (mass) and energy (electronvolt).

Because of this, time is measured in the unusual units of :math:`\text{\AA} \
\sqrt{\text{amu} \ \text{eV}^{-1}}`.

Charge is measured in elementary charges, and temperature is specified in kelvin. It is
important to note that in previous versions, different functions in ASE used different
temperature units. Therefore, ASE functions that took a parameter temperature now take a
parameter temperature_K, which should be used instead.
"""
