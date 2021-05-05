#########################
Writing HDF5 Trajectories
#########################

Writing *narupatools* HDF5 trajectories to file is simple.

.. code-block:: python

   from narupatools.ase.openmm import ASEOpenMMDynamics
   from narupatools.frame.hdf5 import

   dynamics = ASEOpenMMDynamics.from_xml_file("neuraminidase.xml")

   # Add the writer
   writer = add_hdf5_writer(dynamics, filename="mytraj.h5", title="My Trajectory")

   dynamics.run(100)

   # Close the writer
   writer.close()

See :func:`~narupatools.frame.hdf5.hdf5.add_hdf5_writer` for more information.