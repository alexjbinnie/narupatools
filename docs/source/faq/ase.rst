Using with ASE
==============

How do I...
-----------

... create a velocity verlet system?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Given an ASE atoms object, you can create dynamics that can run velocity verlet using :obj:`~narupatools.ase.ASEDynamics.create_velocity_verlet`.

.. testsetup:: velocity_verlet

   from ase import Atoms

   atoms = Atoms()

.. testcode:: velocity_verlet

   from narupatools.ase import ASEDynamics

   # Timestep is in Narupa units (picoseconds), not ASE units.
   dynamics = ASEDynamics.create_velocity_verlet(atoms, timestep=0.1)

... create a system using Langevin dynamics?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Given an ASE atoms object, you can create langevin dynamics directly using :obj:`~narupatools.ase.ASEDynamics.create_langevin`.

.. testsetup:: langevin

   from ase import Atoms

   atoms = Atoms()

.. testcode:: langevin

   from narupatools.ase import ASEDynamics

   # All parameters are in Narupa units
   dynamics = ASEDynamics.create_langevin(atoms, timestep=0.1, friction=1e-2, temperature=300)

... run dynamics through ASE?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once you've obtained an ASEDynamics object, go to :ref:`FAQDynamics` to find out what you can do.

... ignore forces and energies?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An ASE atoms object contains information about the system (positions, velocities, masses, etc.). However, it doesn't know how to calculate forces and energy directly. To do this, you must set a calculator for the atoms object.

Narupatools provides a null calculator, which returns 0 for all forces and energies.

.. testsetup:: null_calculator

   from ase import Atoms

   atoms = Atoms()

.. testcode:: null_calculator

   from narupatools.ase.calculators import NullCalculator

   atoms.calc = NullCalculator()