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

"""Shorthand methods for creating various IMD interactions."""

from typing import Any

import numpy as np
from narupa.imd import ParticleInteraction

from .forces import (
    CONSTANT_INTERACTION_TYPE,
    GAUSSIAN_INTERACTION_TYPE,
    SPRING_INTERACTION_TYPE,
)


def constant_interaction(
    *,
    force: np.typing.ArrayLike,
    particles: np.typing.ArrayLike,
    scale: float = 1.0,
    **kwargs: Any,
) -> ParticleInteraction:
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
    return ParticleInteraction(
        interaction_type=CONSTANT_INTERACTION_TYPE,
        particles=np.array(particles, dtype=int),
        scale=scale,
        force=np.array(force, dtype=float),
        **kwargs,
    )


def spring_interaction(
    *,
    particles: np.typing.ArrayLike,
    scale: float = 1.0,
    position: np.typing.ArrayLike,
    **kwargs: Any,
) -> ParticleInteraction:
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
    return ParticleInteraction(
        interaction_type=SPRING_INTERACTION_TYPE,
        particles=np.array(particles, dtype=int),
        scale=scale,
        position=np.array(position, dtype=float),
        **kwargs,
    )


def gaussian_interaction(
    *,
    particles: np.typing.ArrayLike,
    scale: float = 1.0,
    position: np.typing.ArrayLike,
    **kwargs: Any,
) -> ParticleInteraction:
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
    return ParticleInteraction(
        interaction_type=GAUSSIAN_INTERACTION_TYPE,
        particles=np.array(particles, dtype=int),
        scale=scale,
        position=np.array(position, dtype=float),
        **kwargs,
    )
