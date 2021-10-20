import numpy as np
from narupa.trajectory import FrameData
from pythreejs import (
    Mesh,
    SphereGeometry,
    MeshLambertMaterial,
    CylinderGeometry,
    Object3D,
)
from scipy.spatial.transform import Rotation

from narupatools.frame.fields import ParticlePositions, ParticleElements, BondPairs
from narupatools.physics.atomic import vdw_radius
from narupatools.physics.vector import (
    magnitude,
    vector,
    angle,
    cross_product,
    normalized,
)

colors = {1: "#ffffff", 6: "#7c7582", 7: "#16a1f2", 8: "#f22958", 14: "#F0C8A0"}

CYLINDER_GEO = CylinderGeometry(
    radiusTop=0.5,
    radiusBottom=0.5,
    height=1,
    openEnded=True,
)

SPHERE_GEO = SphereGeometry(radius=0.5, widthSegments=12, heightSegments=12)


def create_atom_mesh(*, position: np.ndarray, element: int) -> Mesh:
    radius = 0.5 * vdw_radius(element)
    return Mesh(
        geometry=SPHERE_GEO,
        material=MeshLambertMaterial(color=colors[element]),
        position=tuple(position),
        scale=(radius, radius, radius),
    )


def create_bond_mesh(point1: np.ndarray, point2: np.ndarray, radius: float) -> Mesh:
    bond_direction = point2 - point1
    geo = CYLINDER_GEO

    mesh = Mesh(geometry=geo, material=MeshLambertMaterial(color="white"))

    theta = angle(vector(0, 1, 0), bond_direction)
    axis = normalized(cross_product(vector(0, 1, 0), bond_direction))

    rot = Rotation.from_rotvec(theta * axis)

    mesh.position = tuple(0.5 * (point1 + point2))
    mesh.rotation = (*rot.as_euler("XYZ"), "XYZ")
    mesh.scale = (radius, magnitude(bond_direction), radius)

    return mesh


def ball_and_stick(frame: FrameData) -> Object3D:
    repr = Object3D()
    for position, element in zip(frame[ParticlePositions], frame[ParticleElements]):
        ball = create_atom_mesh(position=position, element=element)
        repr.add(ball)

    if BondPairs in frame:
        for bond in frame[BondPairs]:
            bond = create_bond_mesh(
                frame[ParticlePositions][bond[0]],
                frame[ParticlePositions][bond[1]],
                0.04,
            )
            repr.add(bond)

    return repr
