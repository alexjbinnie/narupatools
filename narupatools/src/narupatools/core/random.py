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

"""Core random number generation."""

import random
import string


def random_float(min: float = 0.0, max: float = 1.0) -> float:
    """
    Generate a uniform random scalar in a range.

    :param min: Minimum value (inclusive).
    :param max: Maximum value (inclusive).
    """
    return random.random() * (max - min) + min


def random_integer(min: int = 0, max: int = 1) -> int:
    """
    Generate a uniform random integer in a range.

    :param min: Minimum value (inclusive).
    :param max: Maximum value (inclusive).
    """
    return random.randint(min, max)


def random_letter() -> str:
    """Generate a random letter from the current locale."""
    return random.choice(string.ascii_letters)


def random_word(min_length: int = 1, max_length: int = 10) -> str:
    """
    Generate a random word formed of random letters in a mixture of cases.

    :param min_length: Minimum number of characters.
    :param max_length: Maximum number of characters.
    """
    length = random_integer(min_length, max_length)
    return "".join([random_letter() for _ in range(length)])
