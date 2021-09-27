"""Dynamic NGLWidget that shows dynamics."""

from typing import Any

from infinite_sets import everything
from nglview import NGLWidget

from narupatools.core.dynamics import SimulationDynamics
from narupatools.frame.fields import ParticlePositions
from narupatools.nglview._structure import FrameDataStructure
from narupatools.util.timing import throttle


class _DynamicsWidget:
    def __init__(self, dynamics: SimulationDynamics):
        self.widget = NGLWidget()
        self.dynamics = dynamics
        frame = dynamics.get_frame(fields=everything())
        structure = FrameDataStructure(frame)
        self.frame_component = self.widget.add_structure(structure)
        dynamics.on_post_step.add_callback(self._on_dynamics_step)

    def _on_dynamics_step(self, **kwargs: Any) -> None:
        self.refresh()

    @throttle(0.05)
    def refresh(self) -> None:
        frame = self.dynamics.get_frame(fields={ParticlePositions.key})
        positions = ParticlePositions.get(frame)
        self.frame_component.set_coordinates(10.0 * positions)

    def show(self) -> NGLWidget:
        return self.widget


def show_dynamics(dynamics: SimulationDynamics) -> NGLWidget:
    """
    Create an NGLWidget that dynamically shows a narupatools :class:`SimulationDynamics`.

    The widget automatically listens to changes in the underlying dynamics and updates accordingly.

    :param dynamics: Dynamics to listen to.
    :return: NGLWidget dynamically showing dynamics.
    """
    return _DynamicsWidget(dynamics).show()
