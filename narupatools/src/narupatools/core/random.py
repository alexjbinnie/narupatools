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


def random_float(*, minimum: float = 0.0, maximum: float = 1.0) -> float:
    """
    Generate a uniform random scalar in a range.

    :param minimum: Minimum value (inclusive).
    :param maximum: Maximum value (inclusive).
    """
    return random.random() * (maximum - minimum) + minimum


def random_integer(*, minimum: int = 0, maximum: int = 1) -> int:
    """
    Generate a uniform random integer in a range.

    :param minimum: Minimum value (inclusive).
    :param maximum: Maximum value (inclusive).
    """
    return random.randint(minimum, maximum)


def random_letter() -> str:
    """Generate a random letter from the current locale."""
    return random.choice(string.ascii_letters)


def random_word(*, minimum_length: int = 1, maximum_length: int = 10) -> str:
    """
    Generate a random word formed of random letters in a mixture of cases.

    :param minimum_length: Minimum number of characters.
    :param maximum_length: Maximum number of characters.
    """
    length = random_integer(minimum=minimum_length, maximum=maximum_length)
    return "".join([random_letter() for _ in range(length)])
