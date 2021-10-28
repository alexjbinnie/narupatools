from __future__ import annotations

import json
from dataclasses import dataclass
from typing import List, Optional

import numpy as np
from infinite_sets import InfiniteSet
from MDAnalysis.topology.tables import SYMB2Z
from narupa.trajectory import FrameData

from narupatools.frame import (
    BondCount,
    BondPairs,
    ChainCount,
    ChainNames,
    FrameSource,
    ParticleCount,
    ParticleElements,
    ParticleNames,
    ParticleResidues,
    ResidueChains,
    ResidueCount,
    ResidueNames,
)
from MDAnalysis.topology.tables import Z2SYMB, masses

@dataclass
class HDF5Atom:
    """Atom in a HDF5 topology."""

    index: int
    residue: HDF5Residue
    name: Optional[str] = None
    element: Optional[str] = None
    mass: Optional[float] = None

    @property
    def atomic_number(self) -> Optional[int]:
        """Atomic number of the atom."""
        if self.element is not None:
            return SYMB2Z[self.element]
        return None


@dataclass
class HDF5Residue:
    """Residue in a HDF5 topology."""

    index: int
    chain: HDF5Chain
    name: Optional[str] = None


@dataclass
class HDF5Chain:
    """Chain in a HDF5 topology."""

    index: int
    name: Optional[str] = None


class HDF5Topology(FrameSource):
    """Topology read from a HDF5 trajectory."""

    def __init__(self) -> None:
        self.atoms: List[HDF5Atom] = []
        self.residues: List[HDF5Residue] = []
        self.chains: List[HDF5Chain] = []
        self.bonds: List[List[int]] = []

    def _add_chain(self) -> HDF5Chain:
        chain = HDF5Chain(len(self.chains))
        self.chains.append(chain)
        return chain

    def _add_residue(self, chain: HDF5Chain) -> HDF5Residue:
        residue = HDF5Residue(len(self.residues), chain)
        self.residues.append(residue)
        return residue

    def _add_atom(self, residue: HDF5Residue) -> HDF5Atom:
        atom = HDF5Atom(len(self.atoms), residue)
        self.atoms.append(atom)
        return atom

    @property
    def masses(self) -> np.ndarray:
        return np.array([atom.mass for atom in self.atoms])

    @classmethod
    def from_string(cls, string: str) -> HDF5Topology:
        """Create a HDF5 topology from a JSON string."""
        content = json.loads(string)
        topology = cls()
        if "bonds" in content:
            topology.bonds = content["bonds"]
        if "chains" in content:
            for chain_json in content["chains"]:
                chain = topology._add_chain()
                if "name" in chain_json:
                    chain.name = chain_json["name"]
                if "residues" in chain_json:
                    for residue_json in chain_json["residues"]:
                        residue = topology._add_residue(chain)
                        if "name" in residue_json:
                            residue.name = residue_json["name"]
                        if "atoms" in residue_json:
                            for atom_json in residue_json["atoms"]:
                                atom = topology._add_atom(residue)
                                if "element" in atom_json:
                                    atom.element = atom_json["element"]
                                if "name" in atom_json:
                                    atom.name = atom_json["name"]
                                if "mass" in atom_json:
                                    atom.mass = atom_json["mass"]
                                else:
                                    atom.mass = masses[Z2SYMB[atom.atomic_number]]
        return topology

    def get_frame(self, *, fields: InfiniteSet[str]) -> FrameData:  # noqa: D102
        frame = FrameData()
        if ParticleNames in fields:
            frame[ParticleNames] = [atom.name or "" for atom in self.atoms]
        if ParticleElements in fields:
            ParticleElements.set(
                frame, [atom.atomic_number or 0 for atom in self.atoms]
            )
        if ParticleResidues in fields:
            frame[ParticleResidues] = [atom.residue.index for atom in self.atoms]
        if ParticleCount in fields:
            frame[ParticleCount] = len(self.atoms)
        if ResidueNames in fields:
            frame[ResidueNames] = [residue.name or "" for residue in self.residues]
        if ResidueChains in fields:
            frame[ResidueChains] = [residue.chain.index for residue in self.residues]
        if ResidueCount in fields:
            frame[ResidueCount] = len(self.residues)
        if ChainNames in fields:
            frame[ChainNames] = [chain.name or "" for chain in self.chains]
        if ChainCount in fields:
            frame[ChainCount] = len(self.chains)
        if BondCount in fields:
            frame[BondCount] = len(self.bonds)
        if BondPairs in fields:
            frame[BondPairs] = self.bonds
        return frame

    def __repr__(self) -> str:
        return (
            f"<HDF5Topology {len(self.atoms)} atom(s), {len(self.residues)} "
            f"residues(s), {len(self.chains)} chain(s), {len(self.bonds)} bond(s)>"
        )
