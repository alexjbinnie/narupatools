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

import itertools
import random

import pytest
from ase.geometry import cellpar_to_cell
from infinite_sets import InfiniteSet, everything

from narupatools.frame import FrameData
from narupatools.frame.fields import (
    BondCount,
    BondPairs,
    BoxVectors,
    ChainCount,
    ChainNames,
    FrameKey,
    ParticleCharges,
    ParticleCount,
    ParticleElements,
    ParticleForces,
    ParticleMasses,
    ParticleNames,
    ParticlePositions,
    ParticleResidues,
    ParticleTypes,
    ParticleVelocities,
    PotentialEnergy,
    ResidueChains,
    ResidueCount,
    ResidueIds,
    ResidueNames,
)
from narupatools.physics.random import random_float, random_vector
from narupatools.util.random import random_integer, random_word


@pytest.fixture(params=range(20))
def seed(request):
    random.seed(request.param)
    return request.param


@pytest.fixture
def particle_count(seed):
    return random_integer(minimum=1, maximum=100)


@pytest.fixture
def residue_count(seed, particle_count):
    return 1 + particle_count // 5


@pytest.fixture
def chain_count(seed, residue_count):
    return 1 + residue_count // 5


@pytest.fixture
def particle_positions(seed, particle_count):
    return [random_vector(max_magnitude=100.0) for _ in range(particle_count)]


@pytest.fixture
def particle_velocities(seed, particle_count):
    return [random_vector(max_magnitude=100.0) for _ in range(particle_count)]


@pytest.fixture
def particle_forces(seed, particle_count):
    return [random_vector(max_magnitude=100.0) for _ in range(particle_count)]


@pytest.fixture
def particle_charges(seed, particle_count):
    return [random_float(minimum=-5.0, maximum=5.0) for _ in range(particle_count)]


@pytest.fixture
def particle_masses(seed, particle_count):
    return [random_float(minimum=-5.0, maximum=5.0) for _ in range(particle_count)]


@pytest.fixture
def particle_names(seed, particle_count):
    return [
        random_word(minimum_length=2, maximum_length=10).capitalize()
        for _ in range(particle_count)
    ]


@pytest.fixture
def particle_elements(seed, particle_count):
    return [random_integer(minimum=1, maximum=92) for _ in range(particle_count)]


@pytest.fixture
def particle_residues(seed, residue_count, particle_count):
    # Ensure each residue index appears at least once
    ids = itertools.chain(
        range(residue_count),
        (
            random_integer(minimum=0, maximum=residue_count - 1)
            for _ in range(particle_count - residue_count)
        ),
    )

    return sorted(ids)


@pytest.fixture
def particle_types(seed, particle_count):
    return [
        random_word(minimum_length=2, maximum_length=10).capitalize()
        for _ in range(particle_count)
    ]


@pytest.fixture
def residue_names(seed, residue_count):
    return [
        random_word(minimum_length=2, maximum_length=10).capitalize()
        for _ in range(residue_count)
    ]


@pytest.fixture
def residue_ids(seed, residue_count):
    return [str(i + 1) for i in range(residue_count)]


@pytest.fixture
def residue_chains(seed, chain_count, residue_count):
    return sorted(
        random_integer(minimum=0, maximum=chain_count - 1) for _ in range(residue_count)
    )


@pytest.fixture
def box_vectors(seed):
    a = random_float(minimum=50.0, maximum=100.0)
    b = random_float(minimum=50.0, maximum=100.0)
    c = random_float(minimum=50.0, maximum=100.0)
    alpha = random_float(minimum=90.0, maximum=90.0)
    beta = random_float(minimum=90.0, maximum=90.0)
    gamma = random_float(minimum=90.0, maximum=90.0)
    return cellpar_to_cell([a, b, c, alpha, beta, gamma])


@pytest.fixture
def chain_names(seed, chain_count):
    return [
        random_word(minimum_length=1, maximum_length=2).capitalize()
        for _ in range(chain_count)
    ]


@pytest.fixture
def potential_energy(seed):
    return random_float(minimum=-100.0, maximum=100.0)


@pytest.fixture
def bond_count(seed, particle_count):
    if particle_count < 2:
        return 0
    return particle_count // 2


@pytest.fixture
def bond_pairs(seed, bond_count, particle_count):
    if bond_count == 0:
        return []
    else:
        bonds = []
        while len(bonds) < bond_count:
            bond = (
                random_integer(maximum=particle_count - 1),
                random_integer(maximum=particle_count - 1),
            )
            if bond[0] == bond[1]:
                continue
            if bond[0] > bond[1]:
                bond = (bond[1], bond[0])
            if bond not in bonds:
                bonds.append(bond)
        return sorted(bonds)


@pytest.fixture
def make_frame(  # noqa: C901
    seed,
    particle_count,
    particle_charges,
    particle_positions,
    particle_velocities,
    particle_forces,
    particle_names,
    particle_masses,
    particle_elements,
    particle_residues,
    particle_types,
    residue_names,
    residue_chains,
    box_vectors,
    potential_energy,
    residue_count,
    chain_count,
    bond_count,
    bond_pairs,
    residue_ids,
    chain_names,
):
    def make(fields: InfiniteSet[FrameKey]):
        frame = FrameData()
        if ParticleCount in fields:
            frame[ParticleCount] = particle_count
        if ParticleCharges in fields:
            frame[ParticleCharges] = particle_charges
        if ParticlePositions in fields:
            frame[ParticlePositions] = particle_positions
        if ParticleVelocities in fields:
            frame[ParticleVelocities] = particle_velocities
        if ParticleForces in fields:
            frame[ParticleForces] = particle_forces
        if ParticleNames in fields:
            frame[ParticleNames] = particle_names
        if ParticleMasses in fields:
            frame[ParticleMasses] = particle_masses
        if ParticleElements in fields:
            frame[ParticleElements] = particle_elements
        if ParticleResidues in fields:
            frame[ParticleResidues] = particle_residues
        if ParticleTypes in fields:
            frame[ParticleTypes] = particle_types
        if ResidueNames in fields:
            frame[ResidueNames] = residue_names
        if ResidueChains in fields:
            frame[ResidueChains] = residue_chains
        if ResidueCount in fields:
            frame[ResidueCount] = residue_count
        if ChainCount in fields:
            frame[ChainCount] = chain_count
        if BoxVectors in fields:
            frame[BoxVectors] = box_vectors
        if PotentialEnergy in fields:
            frame[PotentialEnergy] = potential_energy
        if BondCount in fields:
            frame[BondCount] = bond_count
        if BondPairs in fields:
            frame[BondPairs] = bond_pairs
        if ResidueIds in fields:
            frame[ResidueIds] = residue_ids
        if ChainNames in fields:
            frame[ChainNames] = chain_names
        return frame

    return make


@pytest.fixture
def frame(seed, make_frame):
    return make_frame(everything())
