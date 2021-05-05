import itertools
import random

import pytest
from ase.geometry import cellpar_to_cell
from infinite_sets import InfiniteSet, everything

from narupatools.core.random import random_integer, random_word
from narupatools.frame import (
    BondCount,
    BondPairs,
    BoxVectors,
    ChainCount,
    ChainNames,
    NarupaFrame,
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
from narupatools.frame.fields import FrameKey
from narupatools.physics.random import random_scalar, random_vector


@pytest.fixture(params=range(20))
def seed(request):
    random.seed(request.param)
    return request.param


@pytest.fixture
def particle_count(seed):
    return random_integer(1, 100)


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
    return [random_scalar(min=-5.0, max=5.0) for _ in range(particle_count)]


@pytest.fixture
def particle_masses(seed, particle_count):
    return [random_scalar(min=-5.0, max=5.0) for _ in range(particle_count)]


@pytest.fixture
def particle_names(seed, particle_count):
    return [
        random_word(min_length=2, max_length=10).capitalize()
        for _ in range(particle_count)
    ]


@pytest.fixture
def particle_elements(seed, particle_count):
    return [random_integer(1, 92) for _ in range(particle_count)]


@pytest.fixture
def particle_residues(seed, residue_count, particle_count):
    # Ensure each residue index appears at least once
    ids = itertools.chain(
        range(residue_count),
        (
            random_integer(0, residue_count - 1)
            for _ in range(particle_count - residue_count)
        ),
    )

    return sorted(ids)


@pytest.fixture
def particle_types(seed, particle_count):
    return [
        random_word(min_length=2, max_length=10).capitalize()
        for _ in range(particle_count)
    ]


@pytest.fixture
def residue_names(seed, residue_count):
    return [
        random_word(min_length=2, max_length=10).capitalize()
        for _ in range(residue_count)
    ]


@pytest.fixture
def residue_ids(seed, residue_count):
    return [str(i + 1) for i in range(residue_count)]


@pytest.fixture
def residue_chains(seed, chain_count, residue_count):
    return sorted(random_integer(0, chain_count - 1) for _ in range(residue_count))


@pytest.fixture
def box_vectors(seed):
    a = random_scalar(min=50.0, max=100.0)
    b = random_scalar(min=50.0, max=100.0)
    c = random_scalar(min=50.0, max=100.0)
    alpha = random_scalar(min=90.0, max=90.0)
    beta = random_scalar(min=90.0, max=90.0)
    gamma = random_scalar(min=90.0, max=90.0)
    return cellpar_to_cell([a, b, c, alpha, beta, gamma])


@pytest.fixture
def chain_names(seed, chain_count):
    return [
        random_word(min_length=1, max_length=2).capitalize() for _ in range(chain_count)
    ]


@pytest.fixture
def potential_energy(seed):
    return random_scalar(min=100.0, max=100.0)


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
                random_integer(max=particle_count - 1),
                random_integer(max=particle_count - 1),
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
        frame = NarupaFrame()
        if ParticleCount in fields:
            ParticleCount.set(frame, particle_count)
        if ParticleCharges in fields:
            ParticleCharges.set(frame, particle_charges)
        if ParticlePositions in fields:
            ParticlePositions.set(frame, particle_positions)
        if ParticleVelocities in fields:
            ParticleVelocities.set(frame, particle_velocities)
        if ParticleForces in fields:
            ParticleForces.set(frame, particle_forces)
        if ParticleNames in fields:
            ParticleNames.set(frame, particle_names)
        if ParticleMasses in fields:
            ParticleMasses.set(frame, particle_masses)
        if ParticleElements in fields:
            ParticleElements.set(frame, particle_elements)
        if ParticleResidues in fields:
            ParticleResidues.set(frame, particle_residues)
        if ParticleTypes in fields:
            ParticleTypes.set(frame, particle_types)
        if ResidueNames in fields:
            ResidueNames.set(frame, residue_names)
        if ResidueChains in fields:
            ResidueChains.set(frame, residue_chains)
        if ResidueCount in fields:
            ResidueCount.set(frame, residue_count)
        if ChainCount in fields:
            ChainCount.set(frame, chain_count)
        if BoxVectors in fields:
            BoxVectors.set(frame, box_vectors)
        if PotentialEnergy in fields:
            PotentialEnergy.set(frame, potential_energy)
        if BondCount in fields:
            BondCount.set(frame, bond_count)
        if BondPairs in fields:
            BondPairs.set(frame, bond_pairs)
        if ResidueIds in fields:
            ResidueIds.set(frame, residue_ids)
        if ChainNames in fields:
            ChainNames.set(frame, chain_names)
        return frame

    return make


@pytest.fixture
def frame(seed, make_frame):
    return make_frame(everything())
