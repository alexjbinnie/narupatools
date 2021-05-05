from typing import Sequence

import numpy as np
from MDAnalysis.core.topologyattrs import TopologyAttr

class Topology:
    def __init__(
        self,
        n_atoms: int = ...,
        n_res: int = ...,
        n_seg: int = ...,
        attrs: Sequence[TopologyAttr] = ...,
        residue_segindex: np.ndarray = ...,
        atom_resindex: np.ndarray = ...,
    ): ...
