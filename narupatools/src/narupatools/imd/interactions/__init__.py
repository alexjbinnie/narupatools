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

"""Interactions that can be applied to narupatools/narupa simulations."""

from ._constant import (
    CONSTANT_INTERACTION_TYPE,
    ConstantInteraction,
    constant_interaction,
)
from ._interaction import Interaction
from ._interactiondata import InteractionData
from ._point import (
    GAUSSIAN_INTERACTION_TYPE,
    SPRING_INTERACTION_TYPE,
    PointInteraction,
    PointInteractionData,
    gaussian_interaction,
    spring_interaction,
)
from ._rigidmotion import (
    RIGIDMOTION_INTERACTION_TYPE,
    RigidMotionInteraction,
    RigidMotionInteractionData,
    rigidmotion_interaction,
)

__all__ = [
    "Interaction",
    "InteractionData",
    "ConstantInteraction",
    "PointInteraction",
    "PointInteractionData",
    "constant_interaction",
    "gaussian_interaction",
    "spring_interaction",
    "CONSTANT_INTERACTION_TYPE",
    "GAUSSIAN_INTERACTION_TYPE",
    "SPRING_INTERACTION_TYPE",
    "RIGIDMOTION_INTERACTION_TYPE",
    "RigidMotionInteraction",
    "RigidMotionInteractionData",
    "rigidmotion_interaction",
]
