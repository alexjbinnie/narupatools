Using with LAMMPS
=================

LAMMPS integration in narupatools wraps around the LAMMPS python wrapper.

How do I...
-----------

... create a simulation from a LAMMPS file?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A LAMMPS file can be read in and turned into a :obj:`narupatools.lammps.LAMMPSSimulation` using :obj:`narupatools.lammps.LAMMPSSimulation.from_file`.


... create dynamics from a LAMMPS file?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The method :obj:`narupatools.lammps.LAMMPSDynamics.from_file` can be used directly, without creating a :obj:`narupatools.lammps.LAMMPSSimulation` first. The simulation can be accessed afterwards through :obj:`narupatools.lammps.LAMMPSSimulation.simulation`.

... create a new simulation from scratch?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To create a new LAMMPS simulation from scratch, use the :obj:`narupatools.lammps.LAMMPSSimulation.create_new`. You must provide what unit system to use, which may be one of the options provided `here <https://docs.lammps.org/units.html>`_. If you want to use "lj" units, then instead of specifying "lj" you must pass in a custom :obj:`narupatools.core.UnitSystem`. This is because "lj" is unitless, and narupatools needs to be told how to convert units to and from actual narupatools units.

... handle a specific LAMMPS error or warning?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

LAMMPS normally just throws generic exceptions, or worst prints them to the output without telling Python. Narupatools wraps each call to a LAMMPS command so that any errors and warnings are raised as either :obj:`narupatools.lammps.LAMMPSError` or :obj:`narupatools.lammps.LAMMPSWarning`. Some specific subclasses of these exist in the :obj:`narupatools.lammps.exceptions` module, such as :obj:`UnknownCommandError` when an unknown LAMMPS command is used.

... gather all of a specific atom property?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Gathering collects a value for each atom across all processors and returns them as an array ordered by atom ID. This can be done using the :obj:`~narupatools.lammps.LAMMPSSimulation.gather_atoms` method, providing one of the predefine atom properties found in :obj:`narupatools.lammps.atom_properties`, such as :obj:`~narupatools.lammps.atom_properties.Position` or :obj:`~narupatools.lammps.atom_properties.Velocity`. Note that this **does not** convert it to narupatools units.

... run LAMMPS through ASE?
^^^^^^^^^^^^^^^^^^^^^^^^^^^

LAMMPS can be run through ASE, allowing ASE to perform integration while delegating to LAMMPS to calculate forces and energies. This can be done using the :obj:`narupatools.ase.lammps.LAMMPSCalculator`.

Of use is the :obj:`narupatools.ase.lammps.atoms_from_lammps_simulation` function, which will create an Atoms object from a simulation and add the calculator automatically.