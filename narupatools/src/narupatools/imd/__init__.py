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

"""Module for code relating to interactive molecular dynamics."""

from ._dynamics import InteractiveSimulationDynamics
from ._feature import InteractionFeature
from ._feature_setclear import SetAndClearInteractionFeature
from .interactions import (
    Interaction,
    constant_interaction,
    gaussian_interaction,
    rigidmotion_interaction,
    spring_interaction,
)

__all__ = [
    "InteractiveSimulationDynamics",
    "Interaction",
    "gaussian_interaction",
    "spring_interaction",
    "constant_interaction",
    "rigidmotion_interaction",
    "InteractionFeature",
    "SetAndClearInteractionFeature",
]
