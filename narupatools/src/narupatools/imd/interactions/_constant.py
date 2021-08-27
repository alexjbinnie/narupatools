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

from __future__ import annotations

from typing import Any, Tuple

import numpy as np

from narupatools.core.properties import float_property, numpy_property
from narupatools.physics.force import mass_weighted_forces
from narupatools.physics.rigidbody import center_of_mass
from narupatools.physics.typing import Vector3Array

from ._interaction import Interaction
from ._interactiondata import InteractionData

CONSTANT_INTERACTION_TYPE = "constant"
"""Key identifying the constant interaction type."""


class ConstantInteractionData(InteractionData):
    """Interaction data for a constant force interaction."""

    @numpy_property(dtype=float)
    def force(self) -> None:
        """Constant force to apply, in kilojoules per mole per nanometer."""
        ...

    @float_property
    def scale(self) -> None:
        """Scaling factor to apply to the force."""
        ...


class ConstantInteraction(Interaction[ConstantInteractionData]):
    """
    Interaction that applies constant force.

    The force is not applied to each particle individually. Instead, the force is distributed such
    that the motion of the particles is the same as if the force has been applied to a composite
    particle located at the center of mass, with mass equal to the total mass.
    """

    @numpy_property(dtype=float)
    def force(self) -> np.ndarray:
        """Constant force to apply, in kilojoules per mole per nanometer."""
        ...

    @float_property
    def scale(self) -> None:
        """Scaling factor to apply to the force."""

    def update(self, interaction: ConstantInteractionData) -> None:  # noqa: D102
        super().update(interaction)
        self.force = interaction.force
        self.scale = interaction.scale

    def calculate_forces_and_energy(self) -> Tuple[Vector3Array, float]:  # noqa: D102
        positions = self.dynamics.positions[self.particle_indices]
        masses = self.dynamics.masses[self.particle_indices]
        center = center_of_mass(positions=positions, masses=masses)
        force = self.force * self.scale
        energy = -np.dot(center, force) * self.scale
        return mass_weighted_forces(force=force, masses=masses), energy


InteractionData.register_interaction_type(
    CONSTANT_INTERACTION_TYPE, ConstantInteractionData
)
Interaction.register_interaction_type(CONSTANT_INTERACTION_TYPE, ConstantInteraction)


def constant_interaction(
    *,
    force: np.typing.ArrayLike,
    particles: np.typing.ArrayLike,
    scale: float = 1.0,
    **kwargs: Any,
) -> ConstantInteractionData:
    r"""
    Interaction that applies a constant force to a set of particles.

    This applies a force :math:`F_C` on the center of mass equal to:

    .. math::
        F_C = k F

    where :math:`F` is the force to apply to the center of mass and :math:`k` is a
    scaling factor.

    The force :math:`F_i` experienced by each particle is given by a mass weighting of
    this force:

    .. math::
        F_i = \frac{m_i}{M} F_C

    where :math:`m_i` is the mass of the i-th particle and :math:`M` is the total mass
    of the set of particles the interaction is applied to.

    :param force: Force :math:`F` to apply to the set of particles, in kilojoules per
                  mole per nanometer.
    :param particles: Set of particles to apply the interaction to.
    :param scale: Dimensionless scaling factor :math:`k` to scale the energy and force
                  applied.
    :param kwargs: Arbitrary values to add to the interaction.
    :return: A particle interaction that represents the given interactive force.
    """
    return ConstantInteractionData(
        interaction_type=CONSTANT_INTERACTION_TYPE,
        particles=np.asarray(particles, dtype=int),
        scale=scale,
        force=np.asarray(force, dtype=float),
        **kwargs,
    )
