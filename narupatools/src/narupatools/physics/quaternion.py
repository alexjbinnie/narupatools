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
Quaternions provided by the quaternion python package.

Importing it through this module supresses the warning it produces without numba, as
we don't use those features and hence don't need it.
"""

import warnings

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from quaternion import (
        as_rotation_matrix,
        as_rotation_vector,
        from_rotation_vector,
        as_quat_array,
        quaternion,
    from_vector_part
    )

__all__ = [
    "quaternion",
    "as_rotation_vector",
    "as_rotation_matrix",
    "from_rotation_vector",
    "as_quat_array",
    "from_vector_part"
]
