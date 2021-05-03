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

import numpy as np
import pytest

from narupatools.frame import (
    ParticleForces,
    ParticleVelocities,
    PotentialEnergy,
    SimulationTotalTime,
)


class VillinDynamicsTests:
    particle_count = 8867
    has_real_forces = True
    timestep = 0.005

    def test_dynamics_step_count(self, dynamics):
        assert dynamics.elapsed_steps == 0
        assert dynamics.total_steps == 0
        dynamics.run(steps=10, block=True)
        assert dynamics.elapsed_steps == 10
        assert dynamics.total_steps == 10
        dynamics.reset()
        assert dynamics.elapsed_steps == 0
        assert dynamics.total_steps == 10
        dynamics.run(steps=5, block=True)
        assert dynamics.elapsed_steps == 5
        assert dynamics.total_steps == 15

    def test_dynamics_time(self, dynamics):
        assert dynamics.elapsed_time == pytest.approx(0 * self.timestep)
        assert dynamics.total_time == pytest.approx(0 * self.timestep)
        dynamics.run(steps=10, block=True)
        assert dynamics.elapsed_time == pytest.approx(10 * self.timestep)
        assert dynamics.total_time == pytest.approx(10 * self.timestep)
        dynamics.reset()
        assert dynamics.elapsed_time == pytest.approx(0 * self.timestep)
        assert dynamics.total_time == pytest.approx(10 * self.timestep)
        dynamics.run(steps=5, block=True)
        assert dynamics.elapsed_time == pytest.approx(5 * self.timestep)
        assert dynamics.total_time == pytest.approx(15 * self.timestep)

    def test_dynamics_velocities(self, dynamics):
        dynamics.run(steps=10, block=True)
        frame = dynamics.get_frame((ParticleVelocities.key))
        velocities = ParticleVelocities.get(frame)
        assert len(velocities) == self.particle_count
        assert any(np.linalg.norm(velocity) > 0.001 for velocity in velocities)

    def test_dynamics_forces(self, dynamics):
        dynamics.run(steps=10, block=True)
        frame = dynamics.get_frame((ParticleForces.key))
        forces = ParticleForces.get(frame)
        assert len(forces) == self.particle_count
        if not self.has_real_forces:
            return
        assert any(np.linalg.norm(force) > 0.001 for force in forces)

    def test_potential_energy(self, dynamics):
        dynamics.run(steps=10, block=True)
        frame = dynamics.get_frame((PotentialEnergy.key))
        potential_energy = PotentialEnergy.get(frame)
        if not self.has_real_forces:
            return
        assert abs(potential_energy) > 0.001

    def test_timestamp(self, dynamics):
        frame0 = dynamics.get_frame((SimulationTotalTime.key))
        time0 = SimulationTotalTime.get(frame0)
        assert time0 == 0

        dynamics.run(steps=15, block=True)

        frame1 = dynamics.get_frame((SimulationTotalTime.key))
        time1 = SimulationTotalTime.get(frame1)
        assert dynamics.total_time == pytest.approx(15 * self.timestep)
        assert time1 == pytest.approx(15 * self.timestep)
