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

import pytest
from MDAnalysis import AtomGroup, Universe
from narupa.trajectory import FrameData
from test_classes.converter import NeuraminidaseTestConverter

from narupatools.frame.fields import (
    BondCount,
    ChainCount,
    ParticleCount,
    ParticlePositions,
    ResidueCount,
)
from narupatools.mdanalysis import (
    mdanalysis_atomgroup_to_frame,
    mdanalysis_universe_to_frame,
)


@pytest.fixture(scope="module")
def neuraminidase_universe(neuraminidase_pdb_filename) -> Universe:
    return Universe(neuraminidase_pdb_filename, guess_bonds=True)


@pytest.fixture(scope="module")
def neuraminidase_frame(neuraminidase_universe) -> FrameData:
    return mdanalysis_universe_to_frame(neuraminidase_universe)


@pytest.fixture(scope="module")
def neuraminidase_selection(neuraminidase_universe) -> AtomGroup:
    return neuraminidase_universe.select_atoms("resid 94:97")  # type: ignore


@pytest.fixture(scope="module")
def neuraminidase_selection_frame(neuraminidase_selection) -> FrameData:
    return mdanalysis_atomgroup_to_frame(neuraminidase_selection)


@pytest.mark.converter
class TestMDAnalysisConverter(NeuraminidaseTestConverter):
    bond_count = 3321

    @pytest.fixture(autouse=True, scope="class")
    def frame(self, neuraminidase_frame):
        return neuraminidase_frame

    @pytest.fixture(autouse=True, scope="class")
    def selection_frame(self, neuraminidase_selection_frame):
        return neuraminidase_selection_frame

    def test_import_works(self, neuraminidase_universe):
        assert isinstance(neuraminidase_universe, Universe)

    def test_selection_particle_count(self, neuraminidase_selection, selection_frame):
        assert len(neuraminidase_selection.atoms) == 36
        assert ParticleCount.get(selection_frame) == 36

    def test_selection_particle_positions(self, selection_frame):
        assert len(ParticlePositions.get(selection_frame)) == 36

    def test_selection_residue_count(self, selection_frame):
        assert ResidueCount.get(selection_frame) == 4

    def test_selection_chain_count(self, selection_frame):
        assert ChainCount.get(selection_frame) == 1

    def test_selection_bond_count(self, selection_frame):
        assert BondCount.get(selection_frame) == 37
