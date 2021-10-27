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

"""Dynamics for running OpenMM through ASE."""

from __future__ import annotations

from os import PathLike
from typing import Union, Optional

from ase.md.md import MolecularDynamics
from simtk.openmm.app import Simulation

from narupatools.ase import ASEDynamics
from narupatools.ase.openmm import openmm_simulation_to_ase_molecular_dynamics
from narupatools.openmm import deserialize_simulation


class ASEOpenMMDynamics(ASEDynamics):
    """Dynamics for running OpenMM through ASE."""

    def __init__(self, *, simulation: Simulation, dynamics: MolecularDynamics):
        """
        Create new dynamics that runs in ASE with OpenMM.

        :param simulation: OpenMM simulation that is used in the dynamics.
        :param dynamics: ASE dynamics to run.
        """
        super().__init__(dynamics)
        self._simulation = simulation

    @property
    def simulation(self) -> Simulation:
        """Underlying OpenMM simulation."""
        return self._simulation

    @staticmethod
    def from_xml_file(path: Union[str, bytes, PathLike], /, *, platform: Optional[str] = None) -> ASEOpenMMDynamics:
        """
        Create ASE dynamics using an OpenMM simulation as a calculator.

        :param path: Path to an XML file of a serialized OpenMM simulation.
        :return: A dynamics object that can be broadcast on a Narupa session.
        """
        with open(path) as file:
            return ASEOpenMMDynamics.from_xml_string(file.read(), platform=platform)

    @staticmethod
    def from_xml_string(string: str, /, *, platform: Optional[str] = None) -> ASEOpenMMDynamics:
        """
        Create ASE dynamics using an OpenMM simulation as a calculator.

        :param string: Contents of an XML file of a serialized OpenMM simulation.
        :return: A dynamics object that can be broadcast on a Narupa session.
        """
        simulation = deserialize_simulation(string, platform=platform)
        return ASEOpenMMDynamics.from_simulation(simulation)

    @staticmethod
    def from_simulation(simulation: Simulation, /) -> ASEOpenMMDynamics:
        """
        Create ASE dynamics using the given OpenMM simulation as a calculator.

        :param simulation: OpenMM simulation to wrap.
        :return: Dynamics object that can be broadcast on a Narupa session.
        """
        dynamics = openmm_simulation_to_ase_molecular_dynamics(simulation)
        return ASEOpenMMDynamics(simulation=simulation, dynamics=dynamics)
