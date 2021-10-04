from narupa.trajectory import FrameData

import narupatools.util.properties as properties
from narupatools.state import SharedStateObject
from narupatools.state.typing import Serializable


class ParticleVisualiser(SharedStateObject):
    @properties.string
    def display_name(self) -> str:
        pass

    @properties.string
    def selection(self) -> str:
        pass

    def visualiser(self) -> Serializable:
        pass

    @properties.boolean
    def hide(self) -> bool:
        pass

    @properties.integer
    def layer(self) -> int:
        pass

    @properties.number
    def priority(self) -> float:
        pass

    def frame(self) -> FrameData:
        pass

    @properties.boolean
    def extend(self) -> bool:
        pass
