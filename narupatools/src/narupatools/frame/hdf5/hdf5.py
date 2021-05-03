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

"""Utilities for adding a HDF5 writer to a simulation."""

from typing import Any

from infinite_sets import everything

from narupatools.imd import InteractiveSimulationDynamics
from .writer import HDF5Writer


def add_hdf5_writer(
    *,
    dynamics: InteractiveSimulationDynamics,
    filename: str,
    title: str,
    n_steps: int = 1
) -> HDF5Writer:
    """
    Add an HDF5 Writer to some ASE dynamics.

    :param dynamics: Dynamics to attach this writer to.
    :param filename: Filename to write to.
    :param title: Optional title for this trajectory.
    :param n_steps: Number of steps to take
    :return: HDF5 writer that is listening to the trajectory.
    """
    writer = HDF5Writer(filename=filename, title=title)

    has_logged_initial = False

    def log_initial(**kwargs: Any) -> None:
        for interaction in dynamics.imd.current_interactions.values():
            writer.create_interaction(
                key=interaction.key,
                interaction=interaction.interaction,
                frame_index=dynamics.total_steps,
                potential_energy=interaction.potential_energy,
                forces=interaction.forces,
            )

        nonlocal has_logged_initial

        if has_logged_initial:
            return

        log_step(**kwargs)
        writer.save_topology(dynamics.get_frame(everything()))
        has_logged_initial = True

    def log_step(**kwargs: Any) -> None:
        if dynamics.total_steps % n_steps > 0:
            return
        writer.save_frame(
            coordinates=dynamics.positions,
            velocities=dynamics.velocities,
            forces=dynamics.forces,
            kineticEnergy=dynamics.kinetic_energy,
            potentialEnergy=dynamics.potential_energy,
            time=dynamics.total_time,
        )
        for interaction in dynamics.imd.current_interactions.values():
            writer.save_interaction(
                key=interaction.key,
                interaction=interaction.interaction,
                frame_index=dynamics.total_steps,
                potential_energy=interaction.potential_energy,
                forces=interaction.forces,
            )

    def log_end_interaction(key: str, **kwargs: Any) -> None:
        writer.end_interaction(key=key, frame_index=dynamics.total_steps)

    dynamics.on_pre_step.add_callback(log_initial)
    dynamics.on_post_step.add_callback(log_step)
    dynamics.imd.on_end_interaction.add_callback(log_end_interaction)

    return writer
