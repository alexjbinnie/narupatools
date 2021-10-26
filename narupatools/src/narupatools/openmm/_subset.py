import copy
from typing import Dict

import numpy as np
from simtk.openmm import (
    CustomAngleForce,
    CustomBondForce,
    CustomExternalForce,
    CustomNonbondedForce,
    System,
)
from simtk.openmm.app import Simulation, Topology

from narupatools.frame import select


def _copy_custombondforce(
    force: CustomBondForce, indices: np.ndarray, index_map: Dict[int, int]
):
    new_force = CustomBondForce(force.getEnergyFunction())
    for pi in range(force.getNumGlobalParameters()):
        new_force.addGlobalParameter(
            force.getGlobalParameterName(pi), force.getGlobalParameterDefaultValue(pi)
        )
    for pi in range(force.getNumPerBondParameters()):
        new_force.addPerBondParameter(force.getPerBondParameterName(pi))
    for bi in range(force.getNumBonds()):
        i, j, params = force.getBondParameters(bi)
        if i in indices and j in indices:
            new_force.addBond(index_map[i], index_map[j], params)
    return new_force


def _copy_customangleforce(
    force: CustomAngleForce, indices: np.ndarray, index_map: Dict[int, int]
):
    new_force = CustomAngleForce(force.getEnergyFunction())
    for pi in range(force.getNumGlobalParameters()):
        new_force.addGlobalParameter(
            force.getGlobalParameterName(pi), force.getGlobalParameterDefaultValue(pi)
        )
    for pi in range(force.getNumPerAngleParameters()):
        new_force.addPerAngleParameter(force.getPerAngleParameterName(pi))
    for ai in range(force.getNumAngles()):
        i, j, k, params = force.getAngleParameters(ai)
        if i in indices and j in indices and k in indices:
            new_force.addAngle(index_map[i], index_map[j], index_map[k], params)
    return new_force


def _copy_customexternalforce(
    force: CustomExternalForce, indices: np.ndarray, index_map: Dict[int, int]
):
    new_force = CustomExternalForce(force.getEnergyFunction())
    for pi in range(force.getNumGlobalParameters()):
        new_force.addGlobalParameter(
            force.getGlobalParameterName(pi), force.getGlobalParameterDefaultValue(pi)
        )
    for pi in range(force.getNumPerParticleParameters()):
        new_force.addPerParticleParameter(force.getPerParticleParameterName(pi))
    for pi in range(force.getNumParticles()):
        i, params = force.getParticleParameters(pi)
        if i in indices:
            new_force.addParticle(index_map[i], params)
    return new_force


def _copy_customnonbondedforce(
    force: CustomNonbondedForce, indices: np.ndarray, index_map: Dict[int, int]
):
    new_force = CustomNonbondedForce(force.getEnergyFunction())
    for pi in range(force.getNumGlobalParameters()):
        new_force.addGlobalParameter(
            force.getGlobalParameterName(pi), force.getGlobalParameterDefaultValue(pi)
        )
    for pi in range(force.getNumPerParticleParameters()):
        new_force.addPerParticleParameter(force.getPerParticleParameterName(pi))
    if force.getNumTabulatedFunctions() > 0:
        raise ValueError("Does not support tabulated functions.")
    if force.getNumFunctions() > 0:
        raise ValueError("Does not support tabulated functions.")
    if force.getNumInteractionGroups() > 0:
        raise ValueError("Does not support interaction groups.")
    new_force.setNonbondedMethod(force.getNonbondedMethod())
    new_force.setCutoffDistance(force.getCutoffDistance())
    new_force.setUseSwitchingFunction(force.getUseSwitchingFunction())
    new_force.setSwitchingDistance(force.getSwitchingDistance())
    new_force.setUseLongRangeCorrection(force.getUseLongRangeCorrection())
    for pi in range(force.getNumParticles()):
        params = force.getParticleParameters(pi)
        if pi in indices:
            new_force.addParticle(params)
    for ei in range(force.getNumExclusions()):
        i, j = force.getExclusionParticles(ei)
        if i in indices and j in indices:
            new_force.addExclusion(index_map[i], index_map[j])
    return new_force


def system_subset(system: System, indices: np.ndarray) -> System:
    indices = np.asarray(indices)

    new_system = System()
    index_map = {}

    j = 0
    for i in range(system.getNumParticles()):
        if i in indices:
            mass = system.getParticleMass(i)
            new_system.addParticle(mass)
            index_map[i] = j
            j += 1

    for fi in range(system.getNumForces()):
        force = system.getForce(fi)
        if isinstance(force, CustomBondForce):
            new_force = _copy_custombondforce(force, indices, index_map)
        elif isinstance(force, CustomAngleForce):
            new_force = _copy_customangleforce(force, indices, index_map)
        elif isinstance(force, CustomExternalForce):
            new_force = _copy_customexternalforce(force, indices, index_map)
        elif isinstance(force, CustomNonbondedForce):
            new_force = _copy_customnonbondedforce(force, indices, index_map)
        else:
            raise ValueError(f"Cannot handle force {force}")
        new_system.addForce(new_force)

    return new_system


def topology_subset(topology: Topology, indices: np.ndarray) -> Topology:
    indices = np.asarray(indices)

    atom_map = {}
    residue_map = {}
    chain_map = {}

    j = 0
    for atom in topology.atoms():
        if atom.index in indices:
            atom_map[atom.index] = j
            j += 1

    j = 0
    for residue in topology.residues():
        if set(atom.index for atom in residue.atoms()) & set(indices):
            residue_map[residue.index] = j
            j += 1

    j = 0
    for chain in topology.chains():
        if set(residue.index for residue in chain.residues()) & residue_map.keys():
            chain_map[chain.index] = j
            j += 1

    new_topology = Topology()

    chains = []
    residues = []
    atoms = []

    for chain in topology.chains():
        if chain.index in chain_map:
            chains.append(new_topology.addChain(chain.id))

    for residue in topology.residues():
        if residue.index in residue_map:
            chain = chains[chain_map[residue.chain.index]]
            residues.append(
                new_topology.addResidue(
                    residue.name, chain, residue.id, residue.insertionCode
                )
            )

    for atom in topology.atoms():
        if atom.index in atom_map:
            residue = residues[residue_map[atom.residue.index]]
            atoms.append(
                new_topology.addAtom(atom.name, atom.element, residue, atom.id)
            )

    for bond in topology.bonds():
        if bond.atom1.index in indices and bond.atom2.index in indices:
            new_topology.addBond(
                atoms[atom_map[bond.atom1.index]],
                atoms[atom_map[bond.atom2.index]],
                bond.type,
                bond.order,
            )

    new_topology.setPeriodicBoxVectors(topology.getPeriodicBoxVectors())
    new_topology.setUnitCellDimensions(topology.getUnitCellDimensions())

    return new_topology


def simulation_subset(simulation: Simulation, indices: np.ndarray) -> Simulation:

    if isinstance(indices, str):
        indices = select(simulation, indices)
    else:
        indices = np.array(indices)

    platform = simulation.context.getPlatform()

    new_simulation = Simulation(
        topology=topology_subset(simulation.topology, indices),
        system=system_subset(simulation.system, indices),
        integrator=copy.deepcopy(simulation.integrator),
        platform=platform,
    )

    state = simulation.context.getState(getPositions=True, getVelocities=True)
    new_simulation.context.setPositions(state.getPositions(asNumpy=True)[indices])
    new_simulation.context.setVelocities(state.getVelocities(asNumpy=True)[indices])

    return new_simulation
