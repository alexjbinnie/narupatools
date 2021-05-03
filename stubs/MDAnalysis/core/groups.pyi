from typing import List

import numpy as np
from MDAnalysis import Universe
from MDAnalysis.core.topologyobjects import TopologyGroup

class Group:
    @property
    def ix(self) -> np.ndarray: ...
    @property
    def indices(self) -> np.ndarray: ...
    def __len__(self) -> int: ...

class ResidueGroup(Group):
    @property
    def segindices(self) -> np.ndarray: ...
    @property
    def resnames(self) -> List[str]: ...
    @property
    def resids(self) -> List[str]: ...

class SegmentGroup(Group):
    @property
    def segids(self) -> np.ndarray: ...

class AtomGroup(Group):
    @property
    def universe(self) -> Universe: ...
    @property
    def n_atoms(self) -> int: ...
    @property
    def n_residues(self) -> int: ...
    @property
    def n_segments(self) -> int: ...
    @property
    def atoms(self) -> AtomGroup: ...
    @property
    def residues(self) -> ResidueGroup: ...
    @property
    def segments(self) -> SegmentGroup: ...
    @property
    def positions(self) -> np.ndarray: ...
    @property
    def resindices(self) -> np.ndarray: ...
    @property
    def names(self) -> np.ndarray: ...
    @property
    def bonds(self) -> TopologyGroup: ...
    @property
    def types(self) -> np.ndarray: ...
    @property
    def masses(self) -> np.ndarray: ...
    @property
    def elements(self) -> np.ndarray: ...
    @property
    def dimensions(self) -> np.ndarray: ...
    @property
    def velocities(self) -> np.ndarray: ...
    @property
    def forces(self) -> np.ndarray: ...
    @property
    def charges(self) -> np.ndarray: ...
