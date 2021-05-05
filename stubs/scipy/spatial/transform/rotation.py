import numpy as np


class Rotation:
    @classmethod
    def from_quat(cls, quat: np.ndarray) -> "Rotation":
        ...

    def as_matrix(self) -> np.ndarray:
        ...
