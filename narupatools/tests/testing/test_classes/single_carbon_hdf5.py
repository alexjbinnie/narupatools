from abc import ABCMeta, abstractmethod
from tempfile import NamedTemporaryFile

import mdtraj
import numpy as np
import pytest
from narupatools.frame.hdf5 import HDF5Trajectory, add_hdf5_writer
from narupatools.imd import InteractiveSimulationDynamics, constant_interaction
from narupatools.physics.vector import vector


class SingleCarbonHDF5Tests(metaclass=ABCMeta):
    """
    Tests for a single carbon test system writing to HDF5

    This system consists of a single carbon-12 atom with a mass of 12 daltons. It is
    initially at (5, 5, 5) nm, within in a (10, 10, 10) nm box. A langevin integrator
    with a timestep of 0.01 ps, temperature of 300 K and friction of 0.01 / ps is used.
    """

    @pytest.fixture
    def hdf5_filename(self):
        file = NamedTemporaryFile(suffix=".h5", delete=False, mode="w")
        yield file.name
        file.close()

    @pytest.fixture()
    @abstractmethod
    def dynamics(self) -> InteractiveSimulationDynamics:  # noqa: PT004
        raise NotImplementedError()

    def test_hdf5_writer(self, dynamics, hdf5_filename):
        writer = add_hdf5_writer(
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
        writer = add_hdf5_writer(
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

    def test_hdf5_writer_multiple_interactions(self, dynamics, hdf5_filename):
        writer = add_hdf5_writer(
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
        dynamics.imd.remove_interaction(key1)
        dynamics.run(20)
        dynamics.imd.remove_interaction(key2)
        dynamics.run(20)
        writer.close()

        traj2 = HDF5Trajectory.load_file(hdf5_filename)

        assert len(traj2.interactions) == 2
        interaction1 = traj2.interactions[key1]
        interaction2 = traj2.interactions[key2]
        assert interaction1.duration == pytest.approx(0.4)
        assert interaction2.duration == pytest.approx(0.4)
        assert (
            interaction1.calculate_work() + interaction2.calculate_work()
            == pytest.approx(traj2.interactions.calculate_work())
        )
