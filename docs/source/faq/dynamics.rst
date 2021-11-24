********
Dynamics
********

.. _FAQDynamics:

A dynamics object wraps various packages' ways of running dynamics in one universal interface.

The creation of a dynamics object depends on what you wish to run it through. See the FAQ for ASE, OpenMM and LAMMPS for more details.

How do I...
-----------


... run dynamics for a given number of steps?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. testsetup:: dynamics_run_steps

   from ase import Atoms
   from narupatools.ase import ASEDynamics

   dynamics = ASEDynamics.create_velocity_verlet(Atoms(), timestep=0.1)

.. testcode:: dynamics_run_steps

   dynamics.run(steps=25)

... run dynamics indefinitely?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. testsetup:: dynamics_run_block

   from ase import Atoms
   from narupatools.ase import ASEDynamics

   dynamics = ASEDynamics.create_velocity_verlet(Atoms(), timestep=0.1)

.. testcode:: dynamics_run_steps

   # Run forever in the background. Set block=True to run in the current thread.
   dynamics.run(block=False)

.. testcleanup:: dynamics_run_steps

   dynamics.stop()

... throttle the speed of the dynamics?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default, dynamics runs as fast it can. For various reasons, you may want to limit how quickly the simulation runs steps.

For example, because the dynamics is running as fast as possible, anything that requires more computing (such as interactions) will slow the system down. It is strongly advised to set a rate on your dynamics for this reason.

To set the minimum time (in seconds) between simulation steps, use :obj:`~narupatools.core.dynamics.SimulationDynamics.playback_interval`:

.. testsetup:: dynamics_set_speed

   from ase import Atoms
   from narupatools.ase import ASEDynamics

   dynamics = ASEDynamics.create_velocity_verlet(Atoms(), timestep=0.1)

.. testcode:: dynamics_set_speed

   # Set minimum time per step to 0.05 seconds
   dynamics.playback_interval = 0.05

To set it in terms of 'rate' (how many steps per second), use :obj:`~narupatools.core.dynamics.SimulationDynamics.playback_rate`:

.. testcode:: dynamics_set_speed

   # Set the maximum rate of dynamics as 20 steps per second
   dynamics.playback_rate = 20

... do something at the start of every step?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All dynamics objects have callbacks that inform you when they are about to perform a step. By adding a **callback**, we can call a function just before each step is run:

.. testsetup:: prestep_callback

   from ase import Atoms
   from narupatools.ase import ASEDynamics

   dynamics = ASEDynamics.create_velocity_verlet(Atoms(), timestep=0.1)

.. testcode:: prestep_callback

   def callback(**kwargs):
       # Perform any calculations here. It is not recommended to modify the dynamics directly.
       pass

   dynamics.on_pre_step.add_callback(callback)

... do something after each step?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All dynamics also has a :obj:`~narupatools.core.dynamics.SimulationDynamics.on_post_step` which acts similarly to above, but is called after each dynamics step.

... reset the simulation to its original state?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A simulation can be reset using either the :obj:`~narupatools.core.dynamics.SimulationDynamics.reset` or :obj:`~narupatools.core.dynamics.SimulationDynamics.restart` commands:

.. testsetup:: prestep_callback

   from ase import Atoms
   from narupatools.ase import ASEDynamics

   dynamics = ASEDynamics.create_velocity_verlet(Atoms(), timestep=0.1)

.. testcode:: prestep_callback

   dynamics.reset()

   dynamics.restart()

... get how many picoseconds have elapsed in the simulation?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The number of picoseconds that the simulation has run can be obtained using :obj:`~narupatools.core.dynamics.SimulationDynamics.elapsed_time`. This number resets to 0 each time you reset the dynamics. If you want the total simulation time including all resets, use :obj:`~narupatools.core.dynamics.SimulationDynamics.total_time`.

If you want the number of steps, then you can use :obj:`~narupatools.core.dynamics.SimulationDynamics.elapsed_steps` and :obj:`~narupatools.core.dynamics.SimulationDynamics.total_steps` analogously to the time.

... get or set the timestep of the simulation?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The timestep in picoseconds of the simulation can be found using :obj:`~narupatools.core.dynamics.SimulationDynamics.timestep`. Depending on the simulation, you may be able to alter this after the dynamics has been created.

... get properties of the simulation?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Dynamics define a number of properties, such as :obj:`~narupatools.core.dynamics.SimulationDynamics.positions`, :obj:`~narupatools.core.dynamics.SimulationDynamics.velocities`, :obj:`~narupatools.core.dynamics.SimulationDynamics.forces`, :obj:`~narupatools.core.dynamics.SimulationDynamics.potential_energy` and :obj:`~narupatools.core.dynamics.SimulationDynamics.kinetic_energy`. The useful fact about these is that they are always in Narupa units, regardless of if the dynamics is running in ASE, OpenMM or LAMMPS. Again, depending on the simulation you may be able to set some of these values as well (such as positions or velocities).

... redistribute velocities using Maxwell-Boltzmann?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Narupatools provides a useful function for reassigning velocities based on a Maxwell-Boltzmann distribution:

.. testsetup:: maxwell_boltzmann

   from ase import Atoms
   from narupatools.ase import ASEDynamics

   dynamics = ASEDynamics.create_velocity_verlet(Atoms(), timestep=0.1)

.. testcode:: maxwell_boltzmann

   from narupatools.physics.thermodynamics import maxwell_boltzmann_velocities

   dynamics.velocities = maxwell_boltzmann_velocities(masses=dynamics.masses, temperature=300)

Interactive Molecular Dynamics
==============================

Most of the dynamics supported (through ASE, OpenMM or LAMMPS) support interactive molecular dynamics. This allows forces to be applied during a simulation.

These dynamics have a :obj:`~narupatools.imd.InteractiveSimulationDynamics.imd` property which provides access to IMD features.

How do I...
-----------

... get the total work done by interactions?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To get the total work done by a force, you can use :obj:`~narupatools.imd.InteractionFeature.total_work`. This includes the work done by interactions that have now ended, in addition to interactions that have finished.

.. testsetup:: interactive_forces

   from narupatools.ase import ASEDynamics
   from ase import Atoms

   dynamics = ASEDynamics.create_velocity_verlet(Atoms(), timestep=0.1)

.. testcode:: interactive_forces

   dynamics.imd.total_work

... get the current interactive forces?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To get the current interactive forces, use :obj:`~narupatools.imd.InteractionFeature.forces`:

.. testsetup:: interactive_forces

   from narupatools.ase import ASEDynamics
   from ase import Atoms

   dynamics = ASEDynamics.create_velocity_verlet(Atoms(), timestep=0.1)

.. testcode:: interactive_forces

   dynamics.imd.forces

This will be a NumPy array of shape (N, 3), where N is the number of atoms in the system.

Note that when asking for the forces of the dynamics object itself, that will include both the forcefield forces **and** the interaction forces.

... get the current interaction energy?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To get the sum of all the current interaction energies, use :obj:`~narupatools.imd.InteractionFeature.potential_energy`:

.. testsetup:: interactive_forces

   from narupatools.ase import ASEDynamics
   from ase import Atoms

   dynamics = ASEDynamics.create_velocity_verlet(Atoms(), timestep=0.1)

.. testcode:: interactive_forces

   dynamics.imd.potential_energy