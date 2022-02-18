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

import math
from typing import Any, Optional, Union

import numpy as np
from numpy.linalg import LinAlgError, inv

from narupatools.override import override
from narupatools.physics.rigidbody import (
    center_of_mass,
    center_of_mass_velocity,
    moment_of_inertia_tensor,
    spin_angular_momentum,
)
from narupatools.physics.transformation import Rotation, Translation
from narupatools.physics.typing import ScalarArray, Vector3, Vector3Array
from narupatools.physics.vector import (
    cross_product_matrix,
    left_vector_triple_product_matrix,
    vector,
)
from narupatools.util import properties

from ._feedback import InteractionFeedback
from ._interaction import Interaction
from ._parameters import InteractionParameters


class RigidMotionInteractionData(InteractionParameters):
    """Interaction data for a rigid motion interaction."""

    @properties.numpy(dtype=float, shape=(3,))
    def translation(self) -> None:
        """Translation to apply, in nanometers."""
        ...

    @properties.unit_quaternion
    def rotation(self) -> None:
        """Rotation to apply, as a unit quaternion."""
        ...

    @properties.number
    def scale(self) -> None:
        """Scale of the interaction."""
        ...


class RigidMotionInteractionFeedback(InteractionFeedback):
    """Interaction feedback for a rigid motion interaction."""

    @properties.numpy(dtype=float, shape=(3,))
    def accumulated_translation(self) -> None:
        """Total translation target has moved, in nanometers."""
        ...

    @properties.unit_quaternion
    def accumulated_rotation(self) -> None:
        """Total rotation target has rotated, as a unit quaternion."""
        ...


class RigidMotionInteraction(Interaction[RigidMotionInteractionData]):
    """Interaction that applies rotations and translations using damped springs."""

    _feedback_type = RigidMotionInteractionFeedback

    def __init__(self, **kwargs: Any):
        self.translation: Optional[Translation] = None
        self.rotation: Optional[Rotation] = None
        self.scale: float = 1.0
        super().__init__(**kwargs)
        self._accumulated_displacement: Vector3 = vector(0, 0, 0)
        self._accumulated_rotation = Rotation.identity

        self._center_of_mass_dirty = True
        self._center_of_mass_velocity_dirty = True
        self._angular_velocity_dirty = True

    @override(Interaction.update)
    def update(self, interaction: RigidMotionInteractionData) -> None:  # noqa: D102
        super().update(interaction)
        try:
            self.translation = Translation(interaction.translation)
        except AttributeError:
            self.translation = None
        try:
            self.rotation = Rotation(interaction.rotation)
        except AttributeError:
            self.rotation = None
        self.scale = interaction.scale

    def _get_particle_inertias(self) -> Vector3Array:
        try:
            particle_inertia = self.dynamics.moments_of_inertia[self.particle_indices]
        except AttributeError:
            particle_inertia = None

        return particle_inertia  # type: ignore[no-any-return]

    def _calculate_angular_velocity(
        self, *, positions: Vector3Array, velocities: Vector3Array, masses: ScalarArray
    ) -> Vector3Array:
        particle_inertia = self._get_particle_inertias()

        inertia_tensor = moment_of_inertia_tensor(
            positions=positions, masses=masses, origin=self.center_of_mass
        )

        if particle_inertia is not None:
            if len(particle_inertia.shape) == 1:
                inertia_tensor += particle_inertia.sum() * np.identity(3)
            else:
                raise ValueError("Non-symmetric inertia not supported.")

        angular_momentum = spin_angular_momentum(
            masses=masses,
            positions=positions,
            velocities=velocities,
            center_of_mass=self.center_of_mass,
            center_of_mass_velocity=self.center_of_mass_velocity,
        )

        try:
            particle_angular_momenta = self.dynamics.angular_momenta[
                self.particle_indices
            ]
            angular_momentum += particle_angular_momenta.sum(axis=-2)
        except AttributeError:
            particle_angular_momenta = None

        try:
            return inv(inertia_tensor) @ angular_momentum
        except LinAlgError:
            return -2.0 * angular_momentum / np.trace(inertia_tensor)  # type: ignore

    @property
    def angular_velocity(self) -> np.ndarray:
        """Angular velocity of selected particles."""
        if self._angular_velocity_dirty:
            self._angular_velocity = self._calculate_angular_velocity(
                positions=self.dynamics.positions[self.particle_indices],
                masses=self.dynamics.masses[self.particle_indices],
                velocities=self.dynamics.velocities[self.particle_indices],
            )
            self._angular_velocity_dirty = False
        return self._angular_velocity

    @property
    def center_of_mass(self) -> np.ndarray:
        """Center of mass of selected particles."""
        if self._center_of_mass_dirty:
            self._center_of_mass = center_of_mass(
                positions=self.dynamics.positions[self.particle_indices],
                masses=self.dynamics.masses[self.particle_indices],
            )
            self._center_of_mass_dirty = False
        return self._center_of_mass

    @property
    def center_of_mass_velocity(self) -> np.ndarray:
        """Center of mass velocity of selected particles."""
        if self._center_of_mass_velocity_dirty:
            self._center_of_mass_velocity = center_of_mass_velocity(
                velocities=self.dynamics.velocities[self.particle_indices],
                masses=self.dynamics.masses[self.particle_indices],
            )
            self._center_of_mass_velocity_dirty = False
        return self._center_of_mass_velocity

    @override(Interaction.calculate_forces_and_energy)
    def calculate_forces_and_energy(self) -> None:  # noqa: D102
        positions = self.dynamics.positions[self.particle_indices]
        masses = self.dynamics.masses[self.particle_indices]

        k = self.scale
        M = masses.sum()

        if not self.mass_weighted:
            k *= M

        gamma = 2 * math.sqrt(k * M)

        com = self.center_of_mass
        com_vel = self.center_of_mass_velocity
        omega = self.angular_velocity

        rotation_matrix = left_vector_triple_product_matrix(omega, omega)

        if self.rotation is not None:
            desired_rotation = self.rotation @ ~self._accumulated_rotation

            theta = desired_rotation.rotation_vector

            angular_acceleration = (k * theta - gamma * omega) / M

            rotation_matrix += cross_product_matrix(angular_acceleration)

        if self.translation is not None:
            desired_translation = self.translation - self._accumulated_displacement
            translation_vector = (k * desired_translation - gamma * com_vel) / M
        else:
            translation_vector = np.zeros(3)

        self._forces = np.zeros((len(self), 3))
        self._torques = np.zeros((len(self), 3))

        self._forces = masses[:, np.newaxis] * (
            (positions - com) @ rotation_matrix.T + translation_vector
        )

        # for i in range(len(self.particle_indices)):

        #    r_i = positions[i] - com
        #    self._forces[i] += masses[i] * (rotation_matrix @ r_i)

        #    if self.rotation is not None:
        #        self._torques += masses[i] * inertias[i] * angular_acceleration

        #    if self.translation is not None:
        #        self._forces[i] += masses[i] * translation_vector

        self._energy = 0.0

    @override(Interaction.on_post_step)
    def on_post_step(self, timestep: float, **kwargs: Any) -> None:  # noqa: D102
        super().on_post_step(timestep=timestep, **kwargs)

        ang_vel = self.angular_velocity
        center_velocity = self.center_of_mass_velocity

        self._accumulated_rotation = (
            Rotation.from_rotation_vector(ang_vel * timestep)
            @ self._accumulated_rotation
        )

        self._accumulated_displacement += center_velocity * timestep

    def create_feeback(self) -> InteractionFeedback:  # noqa: D102
        feedback: RigidMotionInteractionFeedback = super().create_feeback()  # type: ignore
        feedback.accumulated_rotation = self._accumulated_rotation.versor
        feedback.accumulated_translation = self._accumulated_displacement
        return feedback

    @override(Interaction.mark_positions_dirty)
    def mark_positions_dirty(self) -> None:  # noqa: D102
        super().mark_positions_dirty()
        self._center_of_mass_dirty = True
        self._center_of_mass_velocity_dirty = True
        self._angular_velocity_dirty = True

    @override(Interaction.mark_velocities_dirty)
    def mark_velocities_dirty(self) -> None:  # noqa: D102
        super().mark_velocities_dirty()
        self._center_of_mass_dirty = True
        self._center_of_mass_velocity_dirty = True
        self._angular_velocity_dirty = True


RIGIDMOTION_INTERACTION_TYPE = "rigid_motion"
"""Interaction type constant for rigid motion interactions."""

InteractionParameters.register_interaction_type(
    RIGIDMOTION_INTERACTION_TYPE, RigidMotionInteractionData
)
Interaction.register_interaction_type(
    RIGIDMOTION_INTERACTION_TYPE, RigidMotionInteraction
)


def rigidmotion_interaction(
    *,
    particles: np.typing.ArrayLike,
    scale: float = 1.0,
    translation: Optional[np.typing.ArrayLike] = None,
    rotation: Union[Rotation, Optional[np.typing.ArrayLike]] = None,
    **kwargs: Any,
) -> RigidMotionInteractionData:
    r"""
    Create a new rigid motion interaction.

    A rigid motion interaction applies the force:

    .. math:: `F_i = - \frac{m_i}{M} k (t + \theta \times r_i) - \frac{m_i}{M} \gamma (V + \omega \times r_i) + m_i \omega \times (\omega \times r_i)

    This is a combination of a linear spring, a torsional spring and a centripetal force.

    :param particles: Set of particles to apply the force to.
    :param scale: Scaling factor for the interaction.
    :param translation: Translation to apply to the set of particles.
    :param rotation: Rotation to apply to the set of particles.
    """
    if translation is not None:
        kwargs["translation"] = translation
    if rotation is not None:
        if isinstance(rotation, Rotation):
            rotation = rotation.versor
        kwargs["rotation"] = rotation
    return RigidMotionInteractionData(
        interaction_type=RIGIDMOTION_INTERACTION_TYPE,
        particles=np.array(particles, dtype=int),
        scale=scale,
        **kwargs,
    )
