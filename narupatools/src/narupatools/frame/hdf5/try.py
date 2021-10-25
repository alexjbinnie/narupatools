import numpy as np

from narupatools.frame.hdf5._traj import HDF5Trajectory
from narupatools.ase.openmm import ASEOpenMMDynamics

dynamics = ASEOpenMMDynamics.from_xml_file("../../../../../sandbox/nanotube.xml")

with HDF5Trajectory.record(dynamics) as traj:
    dynamics.run(20)

traj.save_to_file(filename="test.h5", overwrite_existing=True)

traj2 = HDF5Trajectory.load_file(filename="test.h5")

print(traj)

print(traj2)

