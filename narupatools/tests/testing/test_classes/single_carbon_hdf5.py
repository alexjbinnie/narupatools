import os
from abc import ABCMeta, abstractmethod
from tempfile import NamedTemporaryFile

import mdtraj
import numpy as np
import pytest

from narupatools.frame.hdf5 import HDF5Trajectory, record_hdf5
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
        writer = record_hdf5(
            dynamics=dynamics, filename=hdf5_filename, title="Test Trajectory"
        )
        dynamics.run(100)
        writer.close()

        traj = mdtraj.load_hdf5(hdf5_filename)
        assert traj.n_frames == 101
        assert traj.n_atoms == 1
        assert traj.timestep == pytest.approx(0.01)
        assert traj.time[-1] == pytest.approx(dynamics.total_time)
        assert traj.xyz[-1] == pytest.approx(dynamics.positions)

        traj2 = HDF5Trajectory.load_file(hdf5_filename)
        assert traj2.positions[-1] == pytest.approx(dynamics.positions)
        assert traj2.velocities[-1] == pytest.approx(dynamics.velocities)
        assert traj2.forces[-1] == pytest.approx(dynamics.forces)
        assert traj2.kinetic_energies[-1] == pytest.approx(dynamics.kinetic_energy)
        assert traj2.potential_energies[-1] == pytest.approx(dynamics.potential_energy)
        assert traj2.times[-1] == pytest.approx(dynamics.total_time)

        assert len(traj2.interactions) == 0

    def test_hdf5_writer_imd(self, dynamics, hdf5_filename):
        writer = record_hdf5(
            dynamics=dynamics, filename=hdf5_filename, title="Test Trajectory"
        )
        dynamics.run(20)

        key = dynamics.imd.add_interaction(
            constant_interaction(force=vector(50.0, 0.0, 0.0), particles=[0])
        )

        dynamics.run(60)

        dynamics.imd.remove_interaction(key)

        dynamics.run(20)

        writer.close()

        traj = mdtraj.load_hdf5(hdf5_filename)
        assert traj.n_frames == 101
        assert traj.n_atoms == 1
        assert traj.timestep == pytest.approx(0.01)
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
        assert interaction.indices == np.array([0])
        assert interaction.start_time == pytest.approx(0.2)
        assert interaction.end_time == pytest.approx(0.8)
        assert interaction.duration == pytest.approx(0.6)
        assert len(interaction.forces) == 61
        assert interaction.calculate_work() == pytest.approx(
            dynamics.imd.total_work, rel=1e-3
        )
        assert interaction.calculate_power() == pytest.approx(
            dynamics.imd.total_work / 0.6, rel=1e-3
        )

        assert traj2.interactions.forces.shape == (101, 1, 3)

    def test_hdf5_writer_multiple_interactions(self, dynamics, hdf5_filename):
        writer = record_hdf5(
            dynamics=dynamics, filename=hdf5_filename, title="Test Trajectory"
        )

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
        writer.close()

        traj2 = HDF5Trajectory.load_file(hdf5_filename)

        assert len(traj2.interactions) == 2
        interaction1 = traj2.interactions[key1]
        interaction2 = traj2.interactions[key2]
        assert interaction1.duration == pytest.approx(0.4)
        assert interaction2.duration == pytest.approx(0.4)
        assert interaction1.frame_range == range(20, 61)
        assert (
            interaction1.calculate_work() + interaction2.calculate_work()
            == pytest.approx(traj2.interactions.calculate_work())
        )

    def test_hdf5_on_reset(self, dynamics, hdf5_filename):
        writer = record_hdf5(
            dynamics=dynamics, filename=hdf5_filename, title="Test Trajectory"
        )
        dynamics.run(25)

        interactions = {}
        dynamics.imd.add_source(interactions)

        key = "interaction.my_interaction"
        interactions[key] = constant_interaction(
            force=vector(50.0, 0.0, 0.0), particles=[0]
        )

        dynamics.run(25)

        dynamics.reset()

        dynamics.run(25)

        interactions.clear()

        dynamics.run(25)

        writer.close()

        assert os.path.exists(hdf5_filename)
        hdf_file_name, hdf_file_ext = os.path.splitext(hdf5_filename)
        hdf_filename2 = hdf_file_name + "-2" + hdf_file_ext
        assert os.path.exists(hdf_filename2)

        traj1 = HDF5Trajectory.load_file(hdf5_filename)
        traj2 = HDF5Trajectory.load_file(hdf_filename2)

        assert len(traj1.interactions) == 1
        assert len(traj2.interactions) == 1

        assert traj1.interactions[key].frame_range == range(25, 51)
        assert traj2.interactions[key].frame_range == range(26)

        assert traj1.times[0] == pytest.approx(0.0)
        assert traj2.times[0] == pytest.approx(0.0)

    def test_hdf5_file_exists(self, dynamics, hdf5_filename):
        with open(hdf5_filename, "a"):
            pass

        with pytest.raises(FileExistsError):
            _ = record_hdf5(
                dynamics=dynamics, filename=hdf5_filename, title="Test Trajectory"
            )

    def test_hdf5_suffix(self, dynamics, hdf5_filename):
        hdf_file_name, hdf_file_ext = os.path.splitext(hdf5_filename)
        hdf5_filename2 = hdf_file_name + "-2" + hdf_file_ext
        hdf5_filename3 = hdf_file_name + "-3" + hdf_file_ext

        with open(hdf5_filename2, "a"):
            pass

        assert not os.path.exists(hdf5_filename3)

        writer = record_hdf5(
            dynamics=dynamics, filename=hdf5_filename, title="Test Trajectory"
        )

        dynamics.run(10)
        dynamics.reset()
        dynamics.run(10)

        writer.close()

        assert os.path.exists(hdf5_filename3)

    def test_hdf5_load_missing_file(self, dynamics, hdf5_filename):
        with pytest.raises(OSError):  # noqa: PT011
            HDF5Trajectory.load_file(hdf5_filename)
