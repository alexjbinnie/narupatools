#################
Starting a Server
#################

In *narupatools*, the idea of a server is replaced by that of a **Session**. This separates cleanly the dynamics (such as :class:`~narupatools.ase.dynamics.ASEDynamics` or :class:`~narupatools.openmm.dynamics.OpenMMDynamics`) from the actual server that you connect to. This makes it easier to run dynamics independently of if you are actually broadcasting a simulation.


The session represents the server you create and which people connect to. The session initially is 'empty' and not showing anything. You then tell the server to show your dynamics. This tells the server to start sending frames from the dynamics. For example:

.. testcode::

   from narupatools.openmm import OpenMMDynamics
   from narupatools.app import Session

   dynamics = OpenMMDynamics.from_xml_file("./nanotube.xml")

   with Session(dynamics) as session:
       print(f"Session started on port {session.port}")

       # Uncomment this to start an infinite loop while running the server
       # session.start_loop()

.. testoutput::

   ...

.. testcleanup::

   dynamics.stop()

We need to use some kind of loop in the script to prevent it from ending. The actual dynamics and server are all running on other threads, so we simply need to keep the main thread alive. The :meth:`~narupatools.app.session.Session.start_loop` function runs an infinite loop, checking each second to see if any of the background threads have crashed out.

Running in a Notebook
---------------------

When running in a notebook, you won't want to use the `with` construct, and instead want to create and close the session in separate cells:

.. testcode::

   session = Session()

.. testcode::

   session.close()

In a notebook, you won't need the to use the :meth:`~narupatools.app.session.Session.start_loop` function, as having the notebook open keeps the other threads going.

Health Checks
-------------

If one of the background threads (like the MD loop, or the loop that is sending frames to the clients) throws an exception, it won't be printed out. To actually check if this has happened, you should call the :math:`~narupatools.app.session.Session.health_check` function. This is called automatically every second when using the :meth:`~narupatools.app.session.Session.start_loop` function.

If you're using a notebook and the server seems to be playing up, use the health check to make sure there hasn't been any exceptions.