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
from abc import ABCMeta, abstractmethod

import numpy as np
import pytest
from infinite_sets import everything

from narupatools.app import Client, Session
from narupatools.frame.fields import ParticlePositions
from narupatools.imd import InteractiveSimulationDynamics, constant_interaction
from narupatools.physics.vector import dot_product, sqr_magnitude, vector, zero_vector
from narupatools.util.timing import wait_for


@pytest.mark.dynamics
class SingleCarbonSystemTests(metaclass=ABCMeta):
    """
    Tests for the single carbon test system.

    This system consists of a single carbon-12 atom with a mass of 12 daltons. It is
    initially at (5, 5, 5) nm, within in a (10, 10, 10) nm box. A velocity verlet
    integrator with a timestep of 0.01 ps is used.
    """

    @pytest.fixture
    @abstractmethod
    def dynamics(self) -> InteractiveSimulationDynamics:  # noqa: PT004
        raise NotImplementedError

    @pytest.mark.parametrize(
        ("position"), [vector(3, 0, 0), vector(0, -7, 0), vector(0, 0, 10)]
    )
    def test_set_position(self, dynamics, position):
        dynamics.positions = np.array([position])
        assert dynamics.positions[0] == pytest.approx(position)

    @pytest.mark.parametrize(
        ("velocity"), [vector(3, 0, 0), vector(0, -7, 0), vector(0, 0, 10)]
    )
    def test_set_velocity(self, dynamics, velocity):
        dynamics.velocities = np.array([velocity])
        assert dynamics.velocities[0] == pytest.approx(velocity)

    @pytest.mark.parametrize(
        ("position"), [vector(3, 0, 0), vector(0, -7, 0), vector(0, 0, 10)]
    )
    @pytest.mark.parametrize(
        ("velocity"), [vector(3, 0, 0), vector(0, -7, 0), vector(0, 0, 10)]
    )
    def test_velocity_verlet(self, dynamics, position, velocity):
        dynamics.positions = np.array([position])
        dynamics.velocities = np.array([velocity])
        dynamics.run(1)
        assert dynamics.positions[0] == pytest.approx(
            position + 0.01 * velocity, rel=5e-2
        )
        assert dynamics.velocities[0] == pytest.approx(velocity, 1e-2)
        assert dynamics.kinetic_energy == pytest.approx(
            0.5 * 12.0 * sqr_magnitude(velocity)
        )
        dynamics.run(4)
        assert dynamics.positions[0] == pytest.approx(position + 0.05 * velocity)
        assert dynamics.velocities[0] == pytest.approx(velocity)
        assert dynamics.kinetic_energy == pytest.approx(
            0.5 * 12.0 * sqr_magnitude(velocity)
        )

    def test_timestep(self, dynamics):
        assert dynamics.timestep == 0.01

    def test_stationary(self, dynamics):
        assert dynamics.positions[0] == pytest.approx(vector(5, 5, 5))
        assert dynamics.velocities[0] == pytest.approx(zero_vector())
        assert dynamics.forces[0] == pytest.approx(zero_vector())
        assert dynamics.masses[0] == pytest.approx(12.0)
        assert dynamics.potential_energy == pytest.approx(0.0)
        assert dynamics.kinetic_energy == pytest.approx(0.0)
        dynamics.run(100)
        assert dynamics.positions[0] == pytest.approx(vector(5, 5, 5))
        assert dynamics.velocities[0] == pytest.approx(zero_vector())
        assert dynamics.forces[0] == pytest.approx(zero_vector())
        assert dynamics.potential_energy == pytest.approx(0.0)
        assert dynamics.kinetic_energy == pytest.approx(0.0)

    @pytest.mark.parametrize(("force"), [1.0, 10.0])
    def test_imd(self, dynamics, force):
        force = vector(force, 0.0, 0.0)
        dynamics.imd.add_interaction(constant_interaction(particles=[0], force=force))
        dynamics.run(100)
        assert len(dynamics.imd.current_interactions) == 1
        assert dynamics.forces[0] == pytest.approx(force)
        elapsed_time = dynamics.elapsed_time
        mass = dynamics.masses[0]
        position = vector(5, 5, 5) + 0.5 / mass * elapsed_time * elapsed_time * force
        velocity = force * elapsed_time / mass
        work = dot_product(force, position - vector(5, 5, 5))
        assert elapsed_time == pytest.approx(1)
        assert dynamics.positions[0] == pytest.approx(position, rel=5e-2)
        assert dynamics.velocities[0] == pytest.approx(velocity, rel=5e-2)
        assert dynamics.imd.total_work == pytest.approx(work, rel=5e-2)
        assert dynamics.kinetic_energy == pytest.approx(
            0.5 * 12.0 * sqr_magnitude(velocity), rel=5e-2
        )
        assert dynamics.potential_energy == pytest.approx(
            -dot_product(force, position), rel=5e-2
        )

    @pytest.mark.parametrize(
        ("position", "velocity", "force"),
        [(vector(3, 0, 0), vector(0, -7, 0), vector(0, 0, 10))],
    )
    def test_velocity_verlet_force(self, dynamics, position, velocity, force):
        dynamics.positions = np.array([position])
        dynamics.velocities = np.array([velocity])
        dt = 0.01
        acceleration = force / dynamics.masses[0]
        dynamics.imd.add_interaction(constant_interaction(particles=[0], force=force))
        dynamics.run(10)
        time = dt * 10
        assert dynamics.positions[0] == pytest.approx(
            position + time * velocity + 0.5 * time * time * acceleration, abs=1e-2
        )
        assert dynamics.velocities[0] == pytest.approx(
            velocity + time * acceleration, rel=1e-2
        )
        dS = time * velocity + 0.5 * time * time * acceleration
        assert dynamics.imd.total_work == pytest.approx(dot_product(force, dS))

    def test_get_frame(self, dynamics):
        frame = dynamics.get_frame(everything())
        assert ParticlePositions.get(frame) == pytest.approx(
            np.array([vector(5, 5, 5)])
        )

    @pytest.mark.session
    def test_client_session(self, dynamics):
        with Session(port=0, run_discovery=False) as session:
            session.show(dynamics)

            session.health_check()

            with Client.connect_to_session(session) as client:
                client.subscribe_to_frames()
                wait_for(lambda: ParticlePositions.key in client.current_frame)

                frame = client.current_frame
                assert ParticlePositions.get(frame) == pytest.approx(
                    np.array([vector(5, 5, 5)])
                )

    @pytest.mark.session
    def test_client_session_imd(self, dynamics):
        with Session(port=0, run_discovery=False) as session:
            session.show(dynamics)

            session.health_check()

            with Client.connect_to_session(session) as client:
                client.subscribe_to_frames()
                wait_for(lambda: ParticlePositions.key in client.current_frame)

                force = vector(1, 0, 0)
                interaction_id = client.start_interaction(
                    constant_interaction(force=force, particles=[0])
                )

                wait_for(lambda: interaction_id in session.shared_state)

                dynamics.run(10)
                assert len(dynamics.imd.current_interactions) == 1
                t = dynamics.timestep * 10
                mass = dynamics.masses[0]
                position = vector(5, 5, 5) + 0.5 / mass * t * t * force
                assert dynamics.positions[0] == pytest.approx(position, rel=1e-2)

                client.stop_interaction(interaction_id)

                wait_for(lambda: interaction_id not in session.shared_state)

                dynamics.run(10)
                assert len(dynamics.imd.current_interactions) == 0
