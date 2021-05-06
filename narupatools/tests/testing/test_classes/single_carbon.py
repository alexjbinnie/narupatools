import time
from abc import ABCMeta, abstractmethod

import numpy as np
import pytest

from narupatools.app import Client, Session
from narupatools.frame import ParticleElements, ParticlePositions
from narupatools.imd import InteractiveSimulationDynamics, constant_interaction
from narupatools.physics.vector import dot_product, sqr_magnitude, vector, zero_vector


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
        raise NotImplementedError()

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
        assert dynamics.positions[0] == pytest.approx(position + 0.01 * velocity)
        assert dynamics.velocities[0] == pytest.approx(velocity)
        assert dynamics.kinetic_energy == pytest.approx(
            0.5 * 12.0 * sqr_magnitude(velocity)
        )
        dynamics.run(4)
        assert dynamics.positions[0] == pytest.approx(position + 0.05 * velocity)
        assert dynamics.velocities[0] == pytest.approx(velocity)
        assert dynamics.kinetic_energy == pytest.approx(
            0.5 * 12.0 * sqr_magnitude(velocity)
        )

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
        time = dynamics.elapsed_time
        mass = dynamics.masses[0]
        position = vector(5, 5, 5) + 0.5 / mass * time * time * force
        velocity = force * time / mass
        work = dot_product(force, position - vector(5, 5, 5))
        assert time == pytest.approx(1)
        assert dynamics.positions[0] == pytest.approx(position)
        assert dynamics.velocities[0] == pytest.approx(velocity)
        assert dynamics.imd.total_work == pytest.approx(work, rel=1e-3)
        assert dynamics.kinetic_energy == pytest.approx(
            0.5 * 12.0 * sqr_magnitude(velocity), rel=1e-3
        )
        assert dynamics.potential_energy == pytest.approx(-dot_product(force, position))

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
        dynamics.run(1)
        assert dynamics.positions[0] == pytest.approx(
            position + dt * velocity + 0.5 * dt * dt * acceleration
        )
        assert dynamics.velocities[0] == pytest.approx(velocity + dt * acceleration)
        dS = dt * velocity + 0.5 * dt * dt * acceleration
        assert dynamics.imd.total_work == pytest.approx(dot_product(force, dS))

    def test_get_frame(self, dynamics):
        frame = dynamics.get_frame({ParticlePositions.key, ParticleElements.key})
        assert ParticlePositions.get(frame) == pytest.approx(
            np.array([vector(5, 5, 5)])
        )

    @pytest.mark.session
    def test_client_session(self, dynamics):
        with Session(port=0) as session:
            dynamics.run(block=False)
            session.show(dynamics)

            time.sleep(0.5)

            session.health_check()

            with Client.connect_to_session(session) as client:
                client.subscribe_to_frames()
                client.wait_until_first_frame(2)

                frame = client.current_frame
                assert ParticlePositions.get(frame) == pytest.approx(
                    np.array([vector(5, 5, 5)])
                )

            dynamics.stop(wait=True)

    @pytest.mark.session
    def test_client_session_imd(self, dynamics):
        with Session(port=0) as session:
            dynamics.run(block=False)
            session.show(dynamics)

            time.sleep(0.5)

            session.health_check()

            with Client.connect_to_session(session) as client:
                client.subscribe_to_frames()
                client.wait_until_first_frame(2)

                dynamics.stop(wait=True)
                force = vector(1, 0, 0)
                interaction_id = client.start_interaction(
                    constant_interaction(force=force, particles=[0])
                )
                time.sleep(0.5)
                dynamics.run(10)
                assert len(dynamics.imd.current_interactions) == 1
                t = dynamics.timestep * 10
                mass = dynamics.masses[0]
                position = vector(5, 5, 5) + 0.5 / mass * t * t * force
                assert dynamics.positions[0] == pytest.approx(position)

                client.stop_interaction(interaction_id)
                time.sleep(0.5)
                dynamics.run(10)
                assert len(dynamics.imd.current_interactions) == 0

            dynamics.stop(wait=True)
