import numpy as np

from narupatools.frame.hdf5._traj import HDF5Trajectory
from narupatools.ase.openmm import ASEOpenMMDynamics
from narupatools.imd import constant_interaction

dynamics = ASEOpenMMDynamics.from_xml_file("../../../../../sandbox/nanotube.xml")

dynamics.run(50)

with HDF5Trajectory.record(dynamics) as traj:
    dynamics.run(50)

    dynamics.imd.add_interaction(constant_interaction(
        particles=[60, 61, 62, 63, 64],
        force=[100, 0, 0]
    ))

    dynamics.run(50)

    dynamics.imd.clear_interactions()

    dynamics.run(50)

traj.save_to_file(filename="test.h5", overwrite_existing=True)

traj2 = HDF5Trajectory.load_file(filename="test.h5")

print(traj)

print(traj2)

print(traj.potential_energies.shape)
print(traj.interactions.potential_energies.shape)
print(list(traj.interactions.values())[0].potential_energies.shape)
print(list(traj.interactions.values())[0].parameters.force.shape)

