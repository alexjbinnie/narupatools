Servers and Sessions
====================

Running a server allows users to connect, view and interact with simulations. In *narupatools*, a server is called a **Session**.

How do I...
-----------

... start a new session?
^^^^^^^^^^^^^^^^^^^^^^^^

To start a new session, simply create it:

.. testcode:: start_session

   from narupatools.app import Session

   session = Session()

.. testcleanup:: start_session

   session.close()

This is an **empty** session - it is not currently showing anything.

... show dynamics through a session?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Given some dynamics object, we can instruct the session to broadcast it using :obj:`~narupatools.app.Session.show()`

.. testsetup:: show_dynamics

   from narupatools.app import Session
   from narupatools.ase import ASEDynamics
   from ase import Atoms

   atoms = Atoms()
   dynamics = ASEDynamics.create_velocity_verlet(atoms, timestep=0.1)
   session = Session()

.. testcode:: show_dynamics

   session.show(dynamics)

.. testcleanup:: show_dynamics

   dynamics.stop()
   session.close()

If we already have the dynamics and want to start the server showing it automatically, we can pass it as an argument to Session:

.. testsetup:: show_dynamics_constructor

   from narupatools.app import Session
   from narupatools.ase import ASEDynamics
   from ase import Atoms

   atoms = Atoms()
   dynamics = ASEDynamics.create_velocity_verlet(atoms, timestep=0.1)

.. testcode:: show_dynamics_constructor

   session = Session(dynamics)

.. testcleanup:: show_dynamics_constructor

   dynamics.stop()
   session.close()

If the dynamics is not playing,

... show a frame through a session?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Showing Narupa FrameData's and any object that can be converted to one (such as an MDAnalysis Universe, OpenMM Simulation or ASE Atoms object) is achieved in the same way as showing dynamics, either by calling show() or by passing it to the constructor of the Session.

.. testsetup:: show_frame

   from narupatools.app import Session
   from narupa.trajectory import FrameData

   frame = FrameData()
   session = Session()

.. testcode:: show_frame

   session.show(frame)

.. testcleanup:: show_frame

   session.close()

... change what a session is showing?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By calling show() with a different object, you can change what a server is showing.

... play, pause or restart what the session is showing?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If the current target (the object the session is showing) supports playing and restarting (such as dynamics), this can be called either directly or through the session:

.. testsetup:: play_pause_reset

   from narupatools.app import Session
   from narupatools.ase import ASEDynamics
   from ase import Atoms

   atoms = Atoms()
   dynamics = ASEDynamics.create_velocity_verlet(atoms, timestep=0.1)

.. testcode:: play_pause_reset

    session = Session(dynamics)

    # These two do the same action
    session.play()
    dynamics.play()

    # Likewise
    session.pause()
    dynamics.pause()

    session.restart()
    dynamics.restart()

.. testcleanup:: play_pause_reset

   dynamics.stop()
   session.close()

... close a session?
^^^^^^^^^^^^^^^^^^^^

To ensure cleanup, when finished with a Session you should always call :obj:`~narupatools.app.Session.close()` when finished:

.. testsetup:: close

   from narupatools.app import Session

   session = Session()

.. testcode:: close

    session.close()

If running in a script and not a notebook, then using the session as a context manager will call close:

.. testsetup:: contextmanager

   from narupatools.app import Session

.. testcode:: contextmanager

    with Session(...):
        # ...
        pass

    # When the with block is left, close() is automatically called

... find out the address and port of the session?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Especially if you've use automatic port selection, you may want to find out what port the server is running on.

.. testsetup:: port_address

   from narupatools.app import Session

   session = Session()

.. testcode:: port_address

    # Get the name of the Session, as it will be displayed to clients
    session.name

    # Get the port of the Session
    session.port

    # Get the address of the Session
    session.address

.. testcleanup:: port_address

   session.close()

... choose what port to run on?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To choose a port, pass the port argument to the session constructor:

.. testsetup:: port

   from narupatools.app import Session

.. testcode:: port

   session = Session(port=44222)

.. testcleanup:: port

   session.close()

By default, the port 38801 is used.

... use a random port?
^^^^^^^^^^^^^^^^^^^^^^

If you don't care what port is used, you can choose to have a random free port chosen by using port = 0.

... check why my session is broken?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When running in notebooks, one of the background threads may crash. Due to how Python works, it won't automatically print it to the screen.

To check that the session is currently working, you can call health_check(). If everything is fine, nothing will happen. If something broke in the background, it will throw this error so you can see what went wrong.

Shared State
============

The shared state contains things such as interactions and other objects

How do I...
-----------

... access the shared state of a Session
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use :obj:`~narupatools.app.Session.shared_state`.

... access the shared state of the Client
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use :obj:`~narupatools.app.Client.shared_state`.