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

"""Point interaction."""

from typing import Any, Protocol, Tuple

import numpy as np

from narupatools.core.partial import partialclass
from narupatools.core.properties import float_property, numpy_property
from narupatools.imd.interactions._interactiondata import InteractionData
from narupatools.physics.force import (
    gaussian_force_and_energy,
    mass_weighted_forces,
    spring_force_and_energy,
)
from narupatools.physics.rigidbody import center_of_mass
from narupatools.physics.typing import Vector3, Vector3Array

from ._interaction import Interaction

SPRING_INTERACTION_TYPE = "spring"
GAUSSIAN_INTERACTION_TYPE = "gaussian"


class OffsetForceFunction(Protocol):
    """Describes a force which applies at an offset."""

    def __call__(self, *, offset: Vector3) -> Tuple[Vector3Array, float]:  # noqa: D102
        ...


class PointInteractionData(InteractionData):
    """Data for a point interaction."""

    @numpy_property(dtype=float)
    def position(self) -> None:
        """Position at which the interaction is applied."""
        ...

    @float_property
    def scale(self) -> None:
        """Scaling factor applied to the force of this interaction."""
        ...


class PointInteraction(Interaction[PointInteractionData]):
    """
    Interaction based on a force applied to a point.

    This interaction encompasses both the gaussian and spring interactions of Narupa.
    It applies a force to the particles as a whole based on the offset from their
    center of mass to a defined point.
    """

    def __init__(self, *, force_func: OffsetForceFunction, **kwargs: Any):
        super().__init__(**kwargs)
        self._force_func = force_func

    def update(self, interaction: PointInteractionData) -> None:  # noqa: D102
        super().update(interaction)
        self.position = interaction.position
        self.interaction_scale = interaction.scale

    def calculate_forces_and_energy(self) -> None:  # noqa: D102
        positions = self.dynamics.positions[self.particle_indices]
        masses = self.dynamics.masses[self.particle_indices]
        center = center_of_mass(positions=positions, masses=masses)
        force, energy = self._force_func(offset=center - self.position)
        force *= self.interaction_scale
        energy *= self.interaction_scale
        self._forces = mass_weighted_forces(force=force, masses=masses)
        self._energy = energy


InteractionData.register_interaction_type(
    GAUSSIAN_INTERACTION_TYPE, PointInteractionData
)
InteractionData.register_interaction_type(SPRING_INTERACTION_TYPE, PointInteractionData)

Interaction.register_interaction_type(
    GAUSSIAN_INTERACTION_TYPE,
    partialclass(PointInteraction, force_func=gaussian_force_and_energy),
)
Interaction.register_interaction_type(
    SPRING_INTERACTION_TYPE,
    partialclass(PointInteraction, force_func=spring_force_and_energy),
)


def gaussian_interaction(
    *,
    particles: np.typing.ArrayLike,
    scale: float = 1.0,
    position: np.typing.ArrayLike,
    **kwargs: Any,
) -> PointInteractionData:
    r"""
    Interaction that applies a gaussian force to a set of particles.

    This applies a force :math:`F_C` on the center of mass equal to:

    .. math::
        F_C = - k (r_C - r_0) \exp(-\frac{(r_C - r_0)^2}{2}))

    where :math:`r_C` is the center of mass of the particles, `r_0` is the position of
    the interaction and :math:`k` is a scaling factor.

    The force :math:`F_i` experienced by each particle is given by a mass weighting of
    this force:

    .. math::
        F_i = \frac{m_i}{M} F_C

    where :math:`m_i` is the mass of the i-th particle and :math:`M` is the total mass
    of the set of particles the interaction is applied to.

    :param position: Position that the gaussian well is centered, in nanometers.
    :param particles: Set of particles to apply the interaction to.
    :param scale: Dimensionless scaling factor :math:`k` to scale the energy and force
                  applied.
    :param kwargs: Arbitrary arguments for the interaction.
    :return: A particle interaction that represents the given interactive force.
    """
    return PointInteractionData(
        interaction_type=GAUSSIAN_INTERACTION_TYPE,
        particles=np.array(particles, dtype=int),
        scale=scale,
        position=np.array(position, dtype=float),
        **kwargs,
    )


def spring_interaction(
    *,
    particles: np.typing.ArrayLike,
    scale: float = 1.0,
    position: np.typing.ArrayLike,
    **kwargs: Any,
) -> PointInteractionData:
    r"""
    Interaction that applies a spring force to a set of particles.

    This applies a force :math:`F_C` on the center of mass equal to:

    .. math::
        F_C = - k (r_C - r_0)

    where :math:`r_C` is the center of mass of the particles, `r_0` is the position of
    the interaction and :math:`k` is a scaling factor.

    The force :math:`F_i` experienced by each particle is given by a mass weighting of
    this force:

    .. math::
        F_i = \frac{m_i}{M} F_C

    where :math:`m_i` is the mass of the i-th particle and :math:`M` is the total mass
    of the set of particles the interaction is applied to.

    :param position: Position that the spring is anchored to, in nanometers.
    :param particles: Set of particles to apply the interaction to.
    :param scale: Dimensionless scaling factor :math:`k` to scale the energy and force
                  applied.
    :param kwargs: Arbitrary arguments for the interaction.
    :return: A particle interaction that represents the given interactive force.
    """
    return PointInteractionData(
        interaction_type=SPRING_INTERACTION_TYPE,
        particles=np.array(particles, dtype=int),
        scale=scale,
        position=np.array(position, dtype=float),
        **kwargs,
    )
