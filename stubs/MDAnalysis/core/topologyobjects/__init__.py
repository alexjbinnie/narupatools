import numpy as np


class TopologyGroup:
    @property
    def indices(self) -> np.ndarray:
        ...

    def __len__(self) -> int:
        ...
