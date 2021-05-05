#################
Starting a Server
#################

In *narupatools*, the idea of a server is replaced by that of a **Session**. This separates cleanly the dynamics (such as :class:`~narupatools.ase.dynamics.ASEDynamics` or :class:`~narupatools.openmm.dynamics.OpenMMDynamics`) from the actual server that you connect to. This makes it easier to run dynamics independently of if you are actually broadcasting a simulation.


The session represents the server you create and which people connect to. The session initially is 'empty' and not showing anything. You then tell the server to show your dynamics. This tells the server to start sending frames from the dynamics. For example:

.. testcode::

   from narupatools.openmm import OpenMMDynamics
   from narupatools.core import NarupaSession

   dynamics = OpenMMDynamics.from_xml_file("./nanotube.xml")

   with NarupaSession.start() as session:
       print(f"Session started on port {session.port}")

       # Start showing the dynamics on the server. Note it still hasn't been started yet
       session.show(dynamics)

       # Start running the dynamics indefinitely
       dynamics.run()

.. testcleanup::

   dynamics.stop()