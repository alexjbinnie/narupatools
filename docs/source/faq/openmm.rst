Using with OpenMM
=================

OpenMM integration in narupatools allows dynamics to be run using an OpenMM simulation object.

How do I...
-----------

... convert an OpenMM simulation to a Narupa FrameData
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

See the FAQ for converting objects. OpenMM states and simulations can be converted to Narupa frame data.

... create dynamics from an OpenMM simulation?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you have an OpenMM simulation object, then an OpenMM dynamics object can be created using :obj:`~narupatools.openmm.OpenMMDynamics.from_simulation`:

.. testsetup:: from_simulation

   from simtk.openmm import Simulation

   simulation = OpenMMSimulation()

.. testcode:: from_simulation

   from narupatools.openmm import OpenMMDynamics

   dynamics = OpenMMDynamics.from_simulation(simulation)

... create dynamics from a serialized OpenMM simulation?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A previously serialized OpenMM simulation can be opened using :obj:`~narupatools.openmm.OpenMMDynamics.from_xml_file`:

.. testcode:: from_xml_file

   from narupatools.openmm import OpenMMDynamics

   dynamics = OpenMMDynamics.from_xml_file("serialized_simulation.xml")

... run dynamics through ASE?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

OpenMM can be run through ASE instead of use directly. This uses ASE to provide the integrator, and purely uses OpenMM to provide forces and energies.

To do this, use :obj:`narupatools.ase.openmm.ASEOpenMMDynamics` instead of :obj:`narupatools.openmm.OpenMMDynamics`. It can be constructed the same way as described above.

... run a Velocity-Verlet integrator in OpenMM?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

OpenMM by default does not come with a Velocity-Verlet integrator. Narupatools provides one, which is based on the documented example for a custom integrator in the OpenMM documentation.

.. testcode:: velocity_verlet

   from narupatools.openmm import VelocityVerletIntegrator

   integrator = VelocityVerletIntegrator(timestep=0.1)
