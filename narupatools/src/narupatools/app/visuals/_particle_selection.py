import numpy as np
import narupatools.util.properties as properties
from narupatools.state import SharedStateObject


class ParticleSelection(SharedStateObject):
    @properties.numpy(dtype=int, shape=(None,))
    def particle_ids(self) -> np.ndarray:
        ...

    @properties.string
    def display_name(self) -> str:
        ...
