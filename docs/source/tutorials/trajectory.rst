####################
Viewing Trajectories
####################

A **Session** can broadcast both molecular dynamics and pre-recorded trajectories. These trajectories can be read in through several packages, including MDAAnalysis and MDTraj.

.. testcode::

   from MDAnalysis import Universe

   from narupatools.app import Session
   from narupatools.frame import TrajectoryPlayback

   # Import a DL_POLY history file using MDAnalysis
   # Replace this with whichever method you would like to create an MDAnalysis Universe from your file
   trajectory = Universe("HISTORY", topology_format="HISTORY", guess_bonds=True)

   # Create a trajectory playback
   playback = TrajectoryPlayback(trajectory=trajectory)

   # Make the trajectory loop
   playback.looping = True

   # Make the trajectory run at 5 frames per second
   playback.playback_rate = 5

   with Session(playback) as session:
       print(f"Session started on port {session.port}")

       # Uncomment this to start an infinite loop while running the server
       # session.start_loop()

.. testoutput::

   ...

.. testcleanup::

   playback.stop()
