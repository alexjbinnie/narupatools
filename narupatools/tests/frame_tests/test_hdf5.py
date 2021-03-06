from pathlib import Path

import pytest

from narupatools.frame.hdf5 import HDF5Trajectory


@pytest.fixture
def nanotube_traj():
    return HDF5Trajectory.load_file(Path(__file__).parent / "nanotube.h5")


def test_trajectory_len(nanotube_traj):
    assert len(nanotube_traj) == 1258


def test_topology_len(nanotube_traj):
    assert len(nanotube_traj.topology.atoms) == 65
    assert len(nanotube_traj.topology.residues) == 1
    assert len(nanotube_traj.topology.chains) == 1
    assert len(nanotube_traj.topology.bonds) == 88


def test_positions_shape(nanotube_traj):
    assert nanotube_traj.positions.shape == (1258, 65, 3)


def test_velocities_shape(nanotube_traj):
    assert nanotube_traj.velocities.shape == (1258, 65, 3)


def test_forces_shape(nanotube_traj):
    assert nanotube_traj.forces.shape == (1258, 65, 3)


def test_times_shape(nanotube_traj):
    assert nanotube_traj.times.shape == (1258,)


def test_kinetic_energy_shape(nanotube_traj):
    assert nanotube_traj.kinetic_energies.shape == (1258,)


def test_potential_energy_shape(nanotube_traj):
    assert nanotube_traj.potential_energies.shape == (1258,)


def test_interactions_len(nanotube_traj):
    assert len(nanotube_traj.interactions) == 7


def test_interactions_forces_shape(nanotube_traj):
    assert nanotube_traj.interactions.forces.shape == (
        1258,
        65,
        3,
    )


def test_interactions_energy_shape(nanotube_traj):
    assert nanotube_traj.interactions.potential_energies.shape == (1258,)


def test_interactions_per_particle_work_shape(nanotube_traj):
    assert nanotube_traj.interactions.calculate_per_particle_work().shape == (65,)


def test_interactions_cumulative_work_shape(nanotube_traj):
    assert nanotube_traj.interactions.calculate_cumulative_work().shape == (1258,)


def test_interactions_work_is_cumulative(nanotube_traj):
    assert (
        nanotube_traj.interactions.calculate_cumulative_work()[-1]
        == nanotube_traj.interactions.calculate_work()
    )


def test_interactions_work_is_sum_per_particle(nanotube_traj):
    assert (
        nanotube_traj.interactions.calculate_per_particle_work().sum()
        == nanotube_traj.interactions.calculate_work()
    )


def test_interaction_frame_range(nanotube_traj):
    for interaction in nanotube_traj.interactions.values():
        assert interaction.frame_indices == pytest.approx(list(interaction.frame_range))


def test_interaction_position_shape(nanotube_traj):
    for interaction in nanotube_traj.interactions.values():
        size = len(interaction.frame_indices)
        assert interaction.interaction_positions.shape == (size, 3)


def test_interaction_scales_shape(nanotube_traj):
    for interaction in nanotube_traj.interactions.values():
        size = len(interaction.frame_indices)
        assert interaction.interaction_scales.shape == (size,)


def test_interaction_forces_shape(nanotube_traj):
    for interaction in nanotube_traj.interactions.values():
        size = len(interaction.frame_indices)
        atoms = len(interaction.indices)
        assert interaction.forces.shape == (size, atoms, 3)


def test_interaction_energy_shape(nanotube_traj):
    for interaction in nanotube_traj.interactions.values():
        size = len(interaction.frame_indices)
        assert interaction.potential_energies.shape == (size,)


def test_interaction_per_particle_work_shape(nanotube_traj):
    for interaction in nanotube_traj.interactions.values():
        size = len(interaction.indices)
        assert interaction.calculate_per_particle_work().shape == (size,)


def test_interaction_work_is_cumulative(nanotube_traj):
    for interaction in nanotube_traj.interactions.values():
        assert (
            interaction.calculate_cumulative_work()[-1] == interaction.calculate_work()
        )
