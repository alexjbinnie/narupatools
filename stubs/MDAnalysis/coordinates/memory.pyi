from typing import Literal, Optional, Union

import numpy as np

class MemoryReader:
    def __init__(
        self,
        coordinate_array: np.ndarray,
        order: Union[
            Literal["afc"],
            Literal["acf"],
            Literal["caf"],
            Literal["fac"],
            Literal["fca"],
            Literal["cfa"],
        ] = ...,
        dimensions: Optional[np.ndarray] = ...,
        dt: float = ...,
        filename: Optional[str] = ...,
        velocities: Optional[np.ndarray] = ...,
        forces: Optional[np.ndarray] = ...,
    ): ...
