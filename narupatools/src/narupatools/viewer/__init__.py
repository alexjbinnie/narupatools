from narupatools.physics.atomic import vdw_radius, atomic_radius, covalent_radius
from narupatools.frame.fields import ParticlePositions
from narupatools.physics.vector import magnitude

colors = {
    1: '#ffffff',
    6: '#7c7582',
    7: '#16a1f2',
    8: '#f22958'
}

def create_atom_mesh(*, position, element):
    return Mesh(geometry=SphereGeometry(radius=0.2 * vdw_radius(element), widthSegments=12, heightSegments=12),
                material=MeshLambertMaterial(color=colors[element]),
                position=position)


def create_bond_mesh(point1, point2, radius):
    bond_direction = point2 - point1
    geo = CylinderGeometry(radiusTop=radius, radiusBottom=radius, height=magnitude(bond_direction), openEnded=True)

    mesh = Mesh(geometry=geo,
                material=MeshLambertMaterial(color='white')
                )

    mesh.position = tuple(0.5 * (point1 + point2))
    mesh.lookAt(tuple(point2))

    return mesh

spheres = []
for position, element in zip(frame.particle_positions, frame.particle_elements):
    ball = create_atom_mesh(position=position, element=element)
    spheres.append(ball)

for bond in frame.bond_pairs:
    bond = create_bond_mesh(frame[ParticlePositions][bond[0]], frame[ParticlePositions][bond[1]], 0.01)
    spheres.append(bond)

key_light = DirectionalLight(color='white', position=[3, 5, 1], intensity=0.5)

c = PerspectiveCamera(position=[0, 1, 1], up=[0, 1, 0], children=[key_light])

scene = Scene(children=[*spheres, c, AmbientLight(color='#777777')], background=None)

orbit = OrbitControls(controlling=c, maxAzimuthAngle=9999, maxDistance=9999, maxZoom=9999, minAzimuthAngle=-9999)

renderer = Renderer(camera=c,
                    scene=scene,
                    alpha=True,
                    clearOpacity=0,
                    controls=[orbit],width=512, height=512)
display(renderer)