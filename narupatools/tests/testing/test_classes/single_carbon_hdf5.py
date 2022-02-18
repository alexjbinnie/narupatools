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
from tempfile import NamedTemporaryFile

import mdtraj
import numpy as np
import pytest

from narupatools.frame.hdf5 import HDF5Trajectory
from narupatools.imd import InteractiveSimulationDynamics, constant_interaction
from narupatools.physics.vector import vector


@pytest.mark.dynamics
class SingleCarbonHDF5Tests(metaclass=ABCMeta):
    """
    Tests for a single carbon test system writing to HDF5

    This system consists of a single carbon-12 atom with a mass of 12 daltons. It is
    initially at (5, 5, 5) nm, within in a (10, 10, 10) nm box. A langevin integrator
    with a timestep of 0.01 ps, temperature of 300 K and friction of 0.01 / ps is used.
    """

    @pytest.fixture
    def hdf5_filename(self):
        file = NamedTemporaryFile(suffix=".h5", delete=True, mode="w")
        filename = file.name
        file.close()
        return filename

    @pytest.fixture
    @abstractmethod
    def dynamics(self) -> InteractiveSimulationDynamics:  # noqa: PT004
        raise NotImplementedError

    def test_hdf5_writer(self, dynamics, hdf5_filename):
        with HDF5Trajectory.record(dynamics, filename=hdf5_filename) as traj:
            dynamics.run(100)

        assert traj.n_frames == 101
        assert traj.n_atoms == 1
        assert traj.times[1] - traj.times[0] == pytest.approx(0.01)
        assert traj.times[-1] == pytest.approx(dynamics.total_time)
        assert traj.positions[-1] == pytest.approx(dynamics.positions)

        traj.close()

        traj2 = HDF5Trajectory.load_file(hdf5_filename)
        assert traj2.positions[-1] == pytest.approx(dynamics.positions)
        assert traj2.velocities[-1] == pytest.approx(dynamics.velocities)
        assert traj2.forces[-1] == pytest.approx(dynamics.forces)
        assert traj2.kinetic_energies[-1] == pytest.approx(dynamics.kinetic_energy)
        assert traj2.potential_energies[-1] == pytest.approx(dynamics.potential_energy)
        assert traj2.times[-1] == pytest.approx(dynamics.total_time)

        assert len(traj2.interactions) == 0

    def test_hdf5_writer_imd(self, dynamics, hdf5_filename):
        with HDF5Trajectory.record(
            dynamics,
            filename=hdf5_filename,
            title="Test Trajectory",
            close_file_after=True,
        ):
            dynamics.run(20)

            key = dynamics.imd.add_interaction(
                constant_interaction(force=vector(50.0, 0.0, 0.0), particles=[0])
            )

            dynamics.run(60)

            dynamics.imd.remove_interaction(key)

            dynamics.run(20)

        traj = mdtraj.load_hdf5(hdf5_filename)
        assert traj.n_frames == 101
        assert traj.n_atoms == 1
        assert traj.time[1] - traj.time[0] == pytest.approx(0.01)
        assert traj.time[-1] == pytest.approx(dynamics.total_time)
        assert traj.xyz[-1] == pytest.approx(dynamics.positions)

        traj2 = HDF5Trajectory.load_file(hdf5_filename)
        assert traj2.positions[-1] == pytest.approx(dynamics.positions)
        assert traj2.velocities[-1] == pytest.approx(dynamics.velocities)
        assert traj2.forces[-1] == pytest.approx(dynamics.forces)
        assert traj2.kinetic_energies[-1] == pytest.approx(dynamics.kinetic_energy)
        assert traj2.potential_energies[-1] == pytest.approx(dynamics.potential_energy)
        assert traj2.times[-1] == pytest.approx(dynamics.total_time)

        assert len(traj2.interactions) == 1
        interaction = traj2.interactions[key]
        assert interaction.particle_indices == np.array([0])
        assert interaction.start_time == pytest.approx(0.2)
        assert interaction.end_time == pytest.approx(0.8)
        assert interaction.duration == pytest.approx(0.6)
        assert len(interaction.forces) == 61
        assert interaction.calculate_work() == pytest.approx(
            dynamics.imd.total_work, rel=1e-3
        )

        assert traj2.interactions.forces.shape == (101, 1, 3)

    def test_hdf5_writer_multiple_interactions(self, dynamics, hdf5_filename):
        with HDF5Trajectory.record(
            dynamics, filename=hdf5_filename, close_file_after=True
        ):
            dynamics.run(20)
            key1 = dynamics.imd.add_interaction(
                constant_interaction(force=vector(50.0, 0.0, 0.0), particles=[0])
            )
            dynamics.run(20)
            key2 = dynamics.imd.add_interaction(
                constant_interaction(force=vector(0.0, 50.0, 0.0), particles=[0])
            )
            dynamics.run(20)
            assert len(dynamics.imd.current_interactions) == 2
            dynamics.imd.remove_interaction(key1)
            dynamics.run(20)
            assert len(dynamics.imd.current_interactions) == 1
            dynamics.imd.remove_interaction(key2)
            dynamics.run(20)
            assert len(dynamics.imd.current_interactions) == 0

        traj2 = HDF5Trajectory.load_file(hdf5_filename)

        assert len(traj2.interactions) == 2
        interaction1 = traj2.interactions[key1]
        interaction2 = traj2.interactions[key2]
        assert interaction1.duration == pytest.approx(0.4)
        assert interaction2.duration == pytest.approx(0.4)
        assert interaction1.frame_range == range(20, 61)

    def test_hdf5_file_exists(self, dynamics, hdf5_filename):
        with open(hdf5_filename, "a"):
            pass

        with pytest.raises(FileExistsError):
            _ = HDF5Trajectory.new_file(filename=hdf5_filename, title="Test Trajectory")

    def test_hdf5_load_missing_file(self, dynamics, hdf5_filename):
        with pytest.raises(OSError):  # noqa: PT011
            HDF5Trajectory.load_file(hdf5_filename)
