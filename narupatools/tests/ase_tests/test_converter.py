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

import ase
import pytest
from ase import Atoms
from narupa.trajectory import FrameData
from test_classes.converter import NeuraminidaseTestConverter

from narupatools.ase import ase_atoms_to_frame
from narupatools.frame.fields import (
    BoxVectors,
    ParticleCount,
    ParticleElements,
    ParticleNames,
    ParticlePositions,
    ParticleResidues,
    BondPairs,
    BondCount,
)


@pytest.fixture(scope="module")
def neuraminidase_atoms(neuraminidase_pdb_filename) -> Atoms:
    return ase.io.read(neuraminidase_pdb_filename)  # type: ignore


@pytest.fixture(scope="module")
def neuraminidase_frame(neuraminidase_atoms) -> FrameData:
    return ase_atoms_to_frame(neuraminidase_atoms)


@pytest.mark.converter
class TestASEConverter(NeuraminidaseTestConverter):
    fields = (
        ParticleCount.key,
        ParticlePositions.key,
        ParticleElements.key,
        ParticleNames.key,
        ParticleResidues.key,
        BoxVectors.key,
    )

    @pytest.fixture(autouse=True, scope="class")
    def frame(self, neuraminidase_frame):
        return neuraminidase_frame

    def test_import_works(self, neuraminidase_atoms):
        assert isinstance(neuraminidase_atoms, Atoms)
        assert len(neuraminidase_atoms) == 3564
