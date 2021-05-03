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

from typing import Collection

import pytest

from narupatools.frame import (
    BondCount,
    BondPairs,
    BoxVectors,
    ChainCount,
    ParticleCount,
    ParticleElements,
    ParticleForces,
    ParticleMasses,
    ParticleNames,
    ParticlePositions,
    ParticleResidues,
    ParticleVelocities,
    ResidueChains,
    ResidueCount,
    ResidueNames,
)


class NeuraminidaseTestConverter:
    bond_count = 3314

    fields: Collection[str] = (
        ParticleCount.key,
        ParticlePositions.key,
        ParticleElements.key,
        ParticleNames.key,
        ParticleResidues.key,
        ResidueCount.key,
        ResidueNames.key,
        ResidueChains.key,
        ChainCount.key,
        BondCount.key,
        BondPairs.key,
        BoxVectors.key,
    )

    optional_fields: Collection[str] = (
        ParticleMasses.key,
        ParticleVelocities.key,
        ParticleForces.key,
    )

    def skip_if_not_converted(self, key):
        if key not in self.fields:
            pytest.skip(f"Converter does not support {key}")

    def test_has_all_fields(self, frame):
        assert set(self.fields).issubset(
            set(frame.arrays.keys()) | set(frame.values.keys())
        )

    def test_has_no_optional_fields(self, frame):
        assert (
            set(self.optional_fields)
            & (set(frame.arrays.keys()) | set(frame.values.keys()))
            == set()
        )

    def test_particle_count(self, frame):
        self.skip_if_not_converted(ParticleCount.key)
        assert ParticleCount.get(frame) == 3564

    def test_particle_positions(self, frame):
        self.skip_if_not_converted(ParticlePositions.key)
        positions = ParticlePositions.get(frame)
        assert len(positions) == 3564
        assert positions[246] == pytest.approx([0.0586, 0.6523, 6.7175])

    def test_particle_elements(self, frame):
        self.skip_if_not_converted(ParticleElements.key)
        elements = ParticleElements.get(frame)
        assert len(elements) == 3564
        assert elements[538] == 6

    def test_particle_names(self, frame):
        self.skip_if_not_converted(ParticleNames.key)
        names = ParticleNames.get(frame)
        assert len(names) == 3564
        assert names[34] == "CB"

    def test_box_vectors(self, frame):
        self.skip_if_not_converted(BoxVectors.key)
        box_vectors = BoxVectors.get(frame)
        assert len(box_vectors) == 3
        assert box_vectors[0] == pytest.approx([18.095, 0.0, 0.0])
        assert box_vectors[1] == pytest.approx([0.0, 18.095, 0.0])
        assert box_vectors[2] == pytest.approx([0.0, 0.0, 18.095])

    def test_residue_count(self, frame):
        self.skip_if_not_converted(ResidueCount.key)
        assert ResidueCount.get(frame) == 751

    def test_residue_names(self, frame):
        self.skip_if_not_converted(ResidueNames.key)
        residue_names = ResidueNames.get(frame)
        assert len(residue_names) == 751
        assert residue_names[0] == "ARG"
        assert residue_names[750] == "HOH"

    def test_bond_count(self, frame):
        self.skip_if_not_converted(BondCount.key)
        assert BondCount.get(frame) == self.bond_count

    def test_bond_pairs(self, frame):
        self.skip_if_not_converted(BondPairs.key)
        bonds = BondPairs.get(frame)
        assert len(bonds) == self.bond_count
