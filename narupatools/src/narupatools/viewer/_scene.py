from typing import Any

from IPython.core.display import display
from narupa.trajectory import FrameData
from pythreejs import (
    AmbientLight,
    ArrowHelper,
    DirectionalLight,
    Object3D,
    OrbitControls,
    PerspectiveCamera,
    Renderer,
    Scene,
)

from narupatools.frame import convert
from narupatools.frame.fields import ParticlePositions, ParticleVelocities
from narupatools.physics.vector import magnitude, normalized
from narupatools.viewer._representation import ball_and_stick


class ViewerScene:
    def __init__(self) -> None:
        self._key_light = DirectionalLight(
            color="white", position=[3, 5, 1], intensity=0.5
        )

        self._camera = PerspectiveCamera(
            position=[0, 1, 1], up=[0, 1, 0], children=[self._key_light]
        )

        self._ambient = AmbientLight(color="#777777")

        self._scene = Scene(children=[self._camera, self._ambient], background=None)

        self._orbit = OrbitControls(
            controlling=self._camera,
            maxAzimuthAngle=9999,
            maxDistance=9999,
            maxZoom=9999,
            minAzimuthAngle=-9999,
        )

        self._renderer = Renderer(
            camera=self._camera,
            scene=self._scene,
            alpha=True,
            clearOpacity=0,
            controls=[self._orbit],
            width=512,
            height=512,
        )

    def display(self) -> Any:
        return display(self._renderer)

    def add_frame(self, frame: Any, render_func: Any) -> None:
        self._scene.add(render_func(frame))

        if ParticleVelocities in frame:
            velocity_arrows = Object3D()
            for vel, pos in zip(frame[ParticleVelocities], frame[ParticlePositions]):
                arrow = ArrowHelper(
                    dir=tuple(normalized(vel)),
                    origin=tuple(pos),
                    length=magnitude(vel),
                    color="#55aaff",
                )
                velocity_arrows.add(arrow)

            self._scene.add(velocity_arrows)


def show(frame: Any) -> ViewerScene:
    if not isinstance(frame, FrameData):
        frame = convert(frame, FrameData)
    viewer = ViewerScene()
    viewer.add_frame(frame, ball_and_stick)
    return viewer
