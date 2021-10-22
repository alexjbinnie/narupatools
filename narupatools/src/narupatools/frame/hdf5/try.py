from narupatools.frame.hdf5._traj import HDF5Trajectory
from narupatools.ase.openmm import ASEOpenMMDynamics

dynamics = ASEOpenMMDynamics.from_xml_file("../../../../../sandbox/nanotube.xml")

with HDF5Trajectory.record(dynamics) as traj:
    dynamics.run(20)

print(traj)