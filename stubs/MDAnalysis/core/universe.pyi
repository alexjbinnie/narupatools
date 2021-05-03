from typing import Any, Optional, Sequence, Union

from MDAnalysis import AtomGroup
from MDAnalysis.core.topology import Topology
from MDAnalysis.core.topologyobjects import TopologyGroup

class Universe:
    def __init__(
        self,
        file: Union[str, Topology],
        guess_bonds: bool = ...,
        format: Optional[str] = ...,
    ): ...
    def select_atoms(self, sel: str) -> AtomGroup: ...
    @property
    def dimensions(self) -> Sequence[float]: ...
    @property
    def atoms(self) -> AtomGroup: ...
    @property
    def bonds(self) -> TopologyGroup: ...
    @property
    def trajectory(self) -> Any: ...
    @trajectory.setter
    def trajectory(self, value: Any) -> None: ...
