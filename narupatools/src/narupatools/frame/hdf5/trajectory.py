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

"""Reader for the NarupaTools HDF5 format."""

from __future__ import annotations

import contextlib
import re
from collections import Mapping
from typing import Iterator, Optional

import numpy as np
import tables
from tables import File, Group, NoSuchNodeError

from narupatools.physics.typing import IntArray, ScalarArray, Vector3Array


class InteractionsView(Mapping):
    """View of the interactions of a HDF5 trajectory."""

    def __init__(self, source: HDF5Trajectory, interactions: Optional[Group] = None):
        self._source = source
        self._interactions = interactions

    def __len__(self) -> int:
        if self._interactions is None:
            return 0
        return self._interactions._v_nchildren

    def __getitem__(self, item: str) -> InteractionView:
        if self._interactions is None:
            raise KeyError
        return InteractionView(self._source, self._interactions._v_children[item])

    def __iter__(self) -> Iterator[str]:
        if self._interactions is None:
            yield from ()
            return
        for child in self._interactions._f_iter_nodes():
            yield child._v_name

    def calculate_work(self) -> float:
        """Calculate the work done by all interactions, in kilojoules per mole."""
        if self._interactions is None:
            return 0.0
        return sum(interaction.calculate_work() for interaction in self.values())


class InteractionView:
    """View of a specific interaction of a HDF5 trajectory."""

    def __init__(self, source: HDF5Trajectory, interaction: Group):
        self._source = source
        self._interaction = interaction
        self._start_index = interaction._v_attrs["startIndex"]
        self._end_index = interaction._v_attrs["endIndex"]
        self._type = interaction._v_attrs["type"]
        self._indices = np.array(interaction.indices, dtype=int)

    @property
    def type(self) -> str:
        """Type of the interaction."""
        return self._type  # type: ignore

    @property
    def start_frame_index(self) -> int:
        """Index of the start frame."""
        return self._start_index  # type: ignore

    @property
    def end_frame_index(self) -> int:
        """Index of the end frame."""
        return self._end_index  # type: ignore

    @property
    def start_time(self) -> float:
        """Start time of the interaction in picoseconds."""
        return self._source._times[self._start_index]  # type: ignore

    @property
    def end_time(self) -> float:
        """End time of the interaction in picoseconds."""
        return self._source._times[self._end_index]  # type: ignore

    @property
    def duration(self) -> float:
        """Duration of the interaction in picoseconds."""
        return self.end_time - self.start_time

    @property
    def indices(self) -> np.ndarray:
        """Indices of the particles affected by the interaction."""
        return self._indices

    @property
    def forces(self) -> Vector3Array:
        """
        Forces on the particles affected by the interaction.

        The forces are in kilojoules per mole per nanometer.
        """
        return np.array(self._interaction.forces)

    @property
    def potential_energies(self) -> ScalarArray:
        """
        Potential energies of the interaction at each frame.

        The energies are in kilojoules per mole.
        """
        return np.array(self._interaction.potentialEnergy)

    @property
    def interaction_positions(self) -> ScalarArray:
        """Position of the interaction at each frame affected by it, in nanometers."""
        return np.array(self._interaction.position)

    @property
    def interaction_scales(self) -> ScalarArray:
        """Scale of the interaction at each frame affected by it."""
        return np.array(self._interaction.scale)

    @property
    def frame_indices(self) -> IntArray:
        """Indices of the frames affected by the interaction."""
        return np.array(self._interaction.frameIndex)

    def calculate_work(self) -> float:
        r"""
        Calculate the work done by an interaction, in kilojoules per mole.

        The work :math:`W` is given by:

        .. math:: W = \sum_i \int_{t_0}^{t_1} F_i(t) \cdot d x_i(t)

        where the sum over :math:`i` is over the particles affected by the interaction.
        In the case of a discretized trajectory, the trapezoidal rule for evaluating
        integrals is used:

        .. math:: W = \sum_i \sum_{n=n_0}^{n_1} \frac{1}{2} (F_i(t + n \delta t) +
                  F_i(t + (n+1) \delta t) \cdot (x_i(t + (n+1) \delta t) -
                  x_i (t + n \delta t))

        where now the second sum is over adjacent timesteps.
        """
        work = 0.0
        for rel_frame, abs_frame in enumerate(
            range(self._start_index, self._end_index)
        ):
            f0 = self._interaction.forces[rel_frame]
            f1 = self._interaction.forces[rel_frame + 1]
            s0 = self._source._positions[abs_frame]
            s1 = self._source._positions[abs_frame + 1]
            work_this_step = 0.0
            for rel_index, abs_index in enumerate(self.indices):
                # Use trapezoidal rule to calculate single step of integral F.dS
                F = 0.5 * (f0[rel_index] + f1[rel_index])
                dS = s1[abs_index] - s0[abs_index]
                work_this_step += np.dot(F, dS)
            work += work_this_step
        return work

    def calculate_power(self) -> float:
        """
        Calculate the average power put in by this interaction.

        The power is in kilojoules per mole per picosecond.
        """
        work = self.calculate_work()
        return work / self.duration


class HDF5Trajectory:
    """Reader for the NarupaTools HDF5 format."""

    def __init__(self, file: File):
        self._file = file
        conventions = re.split(",| ", self._file.root._v_attrs["conventions"])
        if "Pande" not in conventions:
            raise ValueError(
                "Can't read HDF5 trajectory - Pande not specified in conventions."
            )
        mdtraj_conv_version = self._file.root._v_attrs["conventionVersion"]
        if mdtraj_conv_version != "1.1":
            raise ValueError(
                f"Can't read HDF5 trajectory - MDTraj conventions {mdtraj_conv_version}"
                f" unsupported."
            )
        narupatools_conv_version = self._file.root._v_attrs[
            "narupatoolsConventionVersion"
        ]
        if narupatools_conv_version != "1.0":
            raise ValueError(
                f"Can't read HDF5 trajectory - NarupaTools conventions "
                f"{narupatools_conv_version} unsupported."
            )
        self._positions = self._file.root.coordinates
        with contextlib.suppress(NoSuchNodeError):
            self._velocities = self._file.root.velocities
        with contextlib.suppress(NoSuchNodeError):
            self._forces = self._file.root.forces
        with contextlib.suppress(NoSuchNodeError):
            self._times = self._file.root.time
        with contextlib.suppress(NoSuchNodeError):
            self._kinetic_energies = self._file.root.kineticEnergy
        with contextlib.suppress(NoSuchNodeError):
            self._potential_energies = self._file.root.potentialEnergy

        self._interactions: InteractionsView = InteractionsView(self)
        with contextlib.suppress(NoSuchNodeError):
            self._interactions = InteractionsView(self, self._file.root.interactions)

    @classmethod
    def load_file(cls, filename: str) -> HDF5Trajectory:
        """
        Load a trajectory that conforms to the narupatools HDF5 format.

        :param filename: Filename to load.
        """
        file = tables.open_file(filename, mode="r")
        return cls(file)

    @property
    def interactions(self) -> InteractionsView:
        """Interactions applied during this trajectory."""
        return self._interactions

    @property
    def positions(self) -> Vector3Array:
        """Positions of the particles for each frame, in nanometers."""
        return np.array(self._positions)

    @property
    def velocities(self) -> Vector3Array:
        """Velocities of the particles for each frame, in nanometers per picoseconds."""
        return np.array(self._velocities)

    @property
    def forces(self) -> Vector3Array:
        """Forces on the particles for each frame, in nanometers per picoseconds."""
        return np.array(self._forces)

    @property
    def times(self) -> ScalarArray:
        """Time of each frame, in picoseconds."""
        return np.array(self._times)

    @property
    def kinetic_energies(self) -> ScalarArray:
        """Kinetic energy of each frame, in kilojoules per mole."""
        return np.array(self._kinetic_energies)

    @property
    def potential_energies(self) -> ScalarArray:
        """Potential energy of each frame, in kilojoules per mole."""
        return np.array(self._potential_energies)
