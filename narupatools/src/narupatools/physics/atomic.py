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

"""
Atomic data such as radii.

This data is obtained using Wolfram Mathematica's ElementData[].
"""

import numpy as np

_VDW_RADII = np.array(
    [
        np.nan,  # Virtual
        0.120,  # Hydrogen
        0.140,  # Helium
        0.182,  # Lithium
        0.153,  # Beryllium (Mantina et al. 2013)
        0.192,  # Boron (Mantina et al. 2013)
        0.170,  # Carbon
        0.155,  # Nitrogen
        0.152,  # Oxygen
        0.147,  # Fluorine
        0.154,  # Neon
        0.227,  # Sodium
        0.173,  # Magnesium
        0.184,  # Aluminum (Mantina et al. 2013)
        0.210,  # Silicon
        0.180,  # Phosphorus
        0.180,  # Sulfur
        0.175,  # Chlorine
        0.188,  # Argon
        0.275,  # Potassium
        0.231,  # Calcium (Mantina et al. 2013)
        np.nan,  # Scandium
        np.nan,  # Titanium
        np.nan,  # Vanadium
        np.nan,  # Chromium
        np.nan,  # Manganese
        np.nan,  # Iron
        np.nan,  # Cobalt
        0.163,  # Nickel
        0.140,  # Copper
        0.139,  # Zinc
        0.187,  # Gallium
        0.211,  # Germanium (Mantina et al. 2013)
        0.185,  # Arsenic
        0.190,  # Selenium
        0.185,  # Bromine
        0.202,  # Krypton
        0.303,  # Rubidium (Mantina et al. 2013)
        0.250,  # Strontium (Mantina et al. 2013)
        np.nan,  # Yttrium
        np.nan,  # Zirconium
        np.nan,  # Niobium
        np.nan,  # Molybdenum
        np.nan,  # Technetium
        np.nan,  # Ruthenium
        np.nan,  # Rhodium
        0.163,  # Palladium
        0.172,  # Silver
        0.158,  # Cadmium
        0.193,  # Indium
        0.217,  # Tin
        0.206,  # Antimony (Mantina et al. 2013)
        0.206,  # Tellurium
        0.198,  # Iodine
        0.216,  # Xenon
        0.343,  # Cesium (Mantina et al. 2013)
        0.268,  # Barium (Mantina et al. 2013)
        np.nan,  # Lanthanum
        np.nan,  # Cerium
        np.nan,  # Praseodymium
        np.nan,  # Neodymium
        np.nan,  # Promethium
        np.nan,  # Samarium
        np.nan,  # Europium
        np.nan,  # Gadolinium
        np.nan,  # Terbium
        np.nan,  # Dysprosium
        np.nan,  # Holmium
        np.nan,  # Erbium
        np.nan,  # Thulium
        np.nan,  # Ytterbium
        np.nan,  # Lutetium
        np.nan,  # Hafnium
        np.nan,  # Tantalum
        np.nan,  # Tungsten
        np.nan,  # Rhenium
        np.nan,  # Osmium
        np.nan,  # Iridium
        0.175,  # Platinum
        0.166,  # Gold
        0.155,  # Mercury
        0.196,  # Thallium
        0.202,  # Lead
        0.207,  # Bismuth (Mantina et al. 2013)
        0.197,  # Polonium (Mantina et al. 2013)
        0.202,  # Astatine (Mantina et al. 2013)
        0.220,  # Radon (Mantina et al. 2013)
        0.348,  # Francium (Mantina et al. 2013)
        0.283,  # Radium (Mantina et al. 2013)
        np.nan,  # Actinium
        np.nan,  # Thorium
        np.nan,  # Protactinium
        0.186,  # Uranium
        np.nan,  # Neptunium
        np.nan,  # Plutonium
        np.nan,  # Americium
        np.nan,  # Curium
        np.nan,  # Berkelium
        np.nan,  # Californium
        np.nan,  # Einsteinium
        np.nan,  # Fermium
        np.nan,  # Mendelevium
        np.nan,  # Nobelium
        np.nan,  # Lawrencium
        np.nan,  # Rutherfordium
        np.nan,  # Dubnium
        np.nan,  # Seaborgium
        np.nan,  # Bohrium
        np.nan,  # Hassium
        np.nan,  # Meitnerium
        np.nan,  # Darmstadtium
        np.nan,  # Roentgenium
        np.nan,  # Copernicium
        np.nan,  # Nihonium
        np.nan,  # Flerovium
        np.nan,  # Moscovium
        np.nan,  # Livermorium
        np.nan,  # Tennessine
        np.nan,  # Oganesson
    ]
)

_COVALENT_RADII = np.array(
    [
        np.nan,  # Virtual
        0.031,  # Hydrogen
        0.028,  # Helium
        0.128,  # Lithium
        0.096,  # Beryllium
        0.085,  # Boron
        0.076,  # Carbon
        0.071,  # Nitrogen
        0.066,  # Oxygen
        0.057,  # Fluorine
        0.058,  # Neon
        0.166,  # Sodium
        0.141,  # Magnesium
        0.121,  # Aluminum
        0.111,  # Silicon
        0.107,  # Phosphorus
        0.105,  # Sulfur
        0.102,  # Chlorine
        0.106,  # Argon
        0.203,  # Potassium
        0.176,  # Calcium
        0.170,  # Scandium
        0.160,  # Titanium
        0.153,  # Vanadium
        0.139,  # Chromium
        0.139,  # Manganese
        0.132,  # Iron
        0.126,  # Cobalt
        0.124,  # Nickel
        0.132,  # Copper
        0.122,  # Zinc
        0.122,  # Gallium
        0.120,  # Germanium
        0.119,  # Arsenic
        0.120,  # Selenium
        0.120,  # Bromine
        0.116,  # Krypton
        0.220,  # Rubidium
        0.195,  # Strontium
        0.190,  # Yttrium
        0.175,  # Zirconium
        0.164,  # Niobium
        0.154,  # Molybdenum
        0.147,  # Technetium
        0.146,  # Ruthenium
        0.142,  # Rhodium
        0.139,  # Palladium
        0.145,  # Silver
        0.144,  # Cadmium
        0.142,  # Indium
        0.139,  # Tin
        0.139,  # Antimony
        0.138,  # Tellurium
        0.139,  # Iodine
        0.140,  # Xenon
        0.244,  # Cesium
        0.215,  # Barium
        0.207,  # Lanthanum
        0.204,  # Cerium
        0.203,  # Praseodymium
        0.201,  # Neodymium
        0.199,  # Promethium
        0.198,  # Samarium
        0.198,  # Europium
        0.196,  # Gadolinium
        0.194,  # Terbium
        0.192,  # Dysprosium
        0.192,  # Holmium
        0.189,  # Erbium
        0.190,  # Thulium
        0.187,  # Ytterbium
        0.187,  # Lutetium
        0.175,  # Hafnium
        0.170,  # Tantalum
        0.162,  # Tungsten
        0.151,  # Rhenium
        0.144,  # Osmium
        0.141,  # Iridium
        0.136,  # Platinum
        0.136,  # Gold
        0.132,  # Mercury
        0.145,  # Thallium
        0.146,  # Lead
        0.148,  # Bismuth
        0.140,  # Polonium
        0.150,  # Astatine
        0.150,  # Radon
        0.260,  # Francium
        0.221,  # Radium
        0.215,  # Actinium
        0.206,  # Thorium
        0.200,  # Protactinium
        0.196,  # Uranium
        0.190,  # Neptunium
        0.187,  # Plutonium
        0.180,  # Americium
        0.169,  # Curium
        np.nan,  # Berkelium
        np.nan,  # Californium
        np.nan,  # Einsteinium
        np.nan,  # Fermium
        np.nan,  # Mendelevium
        np.nan,  # Nobelium
        np.nan,  # Lawrencium
        np.nan,  # Rutherfordium
        np.nan,  # Dubnium
        np.nan,  # Seaborgium
        np.nan,  # Bohrium
        np.nan,  # Hassium
        np.nan,  # Meitnerium
        np.nan,  # Darmstadtium
        np.nan,  # Roentgenium
        np.nan,  # Copernicium
        np.nan,  # Nihonium
        np.nan,  # Flerovium
        np.nan,  # Moscovium
        np.nan,  # Livermorium
        np.nan,  # Tennessine
        np.nan,  # Oganesson
    ]
)

_ATOMIC_RADII = np.array(
    [
        np.nan,  # Virtual
        0.053,  # Hydrogen
        0.031,  # Helium
        0.167,  # Lithium
        0.112,  # Beryllium
        0.087,  # Boron
        0.067,  # Carbon
        0.056,  # Nitrogen
        0.048,  # Oxygen
        0.042,  # Fluorine
        0.038,  # Neon
        0.190,  # Sodium
        0.145,  # Magnesium
        0.118,  # Aluminum
        0.111,  # Silicon
        0.098,  # Phosphorus
        0.087,  # Sulfur
        0.079,  # Chlorine
        0.071,  # Argon
        0.243,  # Potassium
        0.194,  # Calcium
        0.184,  # Scandium
        0.176,  # Titanium
        0.171,  # Vanadium
        0.166,  # Chromium
        0.161,  # Manganese
        0.156,  # Iron
        0.152,  # Cobalt
        0.149,  # Nickel
        0.145,  # Copper
        0.142,  # Zinc
        0.136,  # Gallium
        0.125,  # Germanium
        0.114,  # Arsenic
        0.103,  # Selenium
        0.094,  # Bromine
        0.087,  # Krypton
        0.265,  # Rubidium
        0.219,  # Strontium
        0.212,  # Yttrium
        0.206,  # Zirconium
        0.198,  # Niobium
        0.190,  # Molybdenum
        0.183,  # Technetium
        0.178,  # Ruthenium
        0.173,  # Rhodium
        0.169,  # Palladium
        0.165,  # Silver
        0.161,  # Cadmium
        0.156,  # Indium
        0.145,  # Tin
        0.133,  # Antimony
        0.123,  # Tellurium
        0.115,  # Iodine
        0.108,  # Xenon
        0.298,  # Cesium
        0.253,  # Barium
        np.nan,  # Lanthanum
        np.nan,  # Cerium
        0.247,  # Praseodymium
        0.206,  # Neodymium
        0.205,  # Promethium
        0.238,  # Samarium
        0.231,  # Europium
        0.233,  # Gadolinium
        0.225,  # Terbium
        0.228,  # Dysprosium
        0.226,  # Holmium
        0.226,  # Erbium
        0.222,  # Thulium
        0.222,  # Ytterbium
        0.217,  # Lutetium
        0.208,  # Hafnium
        0.20,  # Tantalum
        0.193,  # Tungsten
        0.188,  # Rhenium
        0.185,  # Osmium
        0.180,  # Iridium
        0.177,  # Platinum
        0.174,  # Gold
        0.171,  # Mercury
        0.156,  # Thallium
        0.154,  # Lead
        0.143,  # Bismuth
        0.135,  # Polonium
        0.127,  # Astatine
        0.120,  # Radon
        np.nan,  # Francium
        np.nan,  # Radium
        np.nan,  # Actinium
        np.nan,  # Thorium
        np.nan,  # Protactinium
        np.nan,  # Uranium
        np.nan,  # Neptunium
        np.nan,  # Plutonium
        np.nan,  # Americium
        np.nan,  # Curium
        np.nan,  # Berkelium
        np.nan,  # Californium
        np.nan,  # Einsteinium
        np.nan,  # Fermium
        np.nan,  # Mendelevium
        np.nan,  # Nobelium
        np.nan,  # Lawrencium
        np.nan,  # Rutherfordium
        np.nan,  # Dubnium
        np.nan,  # Seaborgium
        np.nan,  # Bohrium
        np.nan,  # Hassium
        np.nan,  # Meitnerium
        np.nan,  # Darmstadtium
        np.nan,  # Roentgenium
        np.nan,  # Copernicium
        np.nan,  # Nihonium
        np.nan,  # Flerovium
        np.nan,  # Moscovium
        np.nan,  # Livermorium
        np.nan,  # Tennessine
        np.nan,  # Oganesson
    ]
)

def atomic_radius(value: int, /) -> float:
    """Get the atomic radii of one or more atomic elements in nanometers."""
    return _ATOMIC_RADII[value]  # type: ignore


def vdw_radius(value: int, /) -> float:
    """Get the VdW radii of one or more atomic elements in nanometers."""
    return _VDW_RADII[value]  # type: ignore


def covalent_radius(value: int, /) -> float:
    """Get the covalent radii of one or more atomic elements in nanometers."""
    return _COVALENT_RADII[value]  # type: ignore
