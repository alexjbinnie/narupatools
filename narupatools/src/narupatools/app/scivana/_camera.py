import numpy as np

from narupatools.state import SharedStateObject
from narupatools.util import properties


class CameraView(SharedStateObject):
    """Fixed camera view from a 2D application."""

    @properties.auto
    def transformation(self) -> np.ndarray:
        """Transformation of the camera as an affine transformation."""

    @properties.auto
    def display_name(self) -> str:
        """Display name of the camera."""

