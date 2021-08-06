######
LAMMPS
######

LAMMPS is a powerful and extensible molecular dynamics code which can be used for a variety of purposes. For all information about LAMMPS, their extensive and thorough `documentation <https://lammps.sandia.gov/doc/Manual.html>`_ is highly recommended.

This will assume that you have installed LAMMPS correctly in your environment, including the python bindings.

Ways of using LAMMPS
====================

LAMMPS operates on specific input files, normally named with a prefix 'in.' such as 'in.peptide'. For interfacing with python there are two ways of using LAMMPS.

The first is what is used by the narupa-lammps package, and allows small snippets of python code to be executed at various points in the LAMMPS MD loop.

The second is invoking LAMMPS commands from python. Here, python is the outermost layer that will control program execution. *narupatools* provides means to work with LAMMPS in this way.

To use your input files with *narupatools*, the only change required would be to remove running logic, such as 'run 100', that appears after the system setup. In both cases presented below, python will be calling 'run' directly to drive the LAMMPS simulation.

Creating a Simulation
=====================

LAMMPS provides a handy python wrapper around LAMMPS, known as PyLammps. For our purposes, we have an additional wrapper known as a :class:`~narupatools.lammps.simulation.LAMMPSSimulation`. This provides several benefits, such as:

* Taking errors and warnings from LAMMPS (which are printed to the console with the prefixes 'ERROR:' and 'WARNING:' respectively) and converting them to Python warnings and exceptions for a more consistent user experience.
* Automatically adding a 'atom_modify map yes' at the start of the file, to ensure atom IDs exist and hence gathering and scattering per-atom data works correctly.
* Accounting for the unit system of your simulation and exposing it in a consistent way. All methods and properties of :class:`~narupatools.lammps.simulation.LAMMPSSimulation` are in standard Narupa/narupatools units, and converted automatically based on the unit system specified in the LAMMPS input file.
* Tries to provide more detailed errors when things go wrong.

A :class:`~narupatools.lammps.simulation.LAMMPSSimulation` should not be created directly. Instead, either:

* Call :meth:`~narupatools.lammps.simulation.LAMMPSSimulation.from_file`, passing in your input file and your data file. The data file is passed in separately as it is read into an MDAnalysis universe to store topology information.
* Call :meth:`~narupatools.lammps.simulation.LAMMPSSimulation.create_new`. This is less recommended, and is mainly used for creating test systems.

As an ASE Calculator
====================

*narupatools* provides a calculator that uses LAMMPS to calculate properties such as forces and energies, whilst allowing ASE to run the dynamics directly. It achieves this in a similar way to the :class:`~narupatools.openmm.OpenMMCalculator`. First, a LAMMPS simulation must be created. Then an ASE atoms object can be created that will be assigned a LAMMPS calcuator, using the :func:`~narupatools.lammps.atoms_from_lammps_simulation` function.

This :class:`~ase.Atoms` object can now be used with minimization and molecular dynamics routines within LAMMPS. When the LAMMPS calculator is asked for forces or potential energy, it first ensures that the latest positions and velocities are copied from the ASE atoms object to the underlying LAMMPS simulation. It then runs a 'run 0' command, to update the simulation within LAMMPS. Finally, the required forces and energies are extracted from LAMMPS.

For example, to use the peptide example provided with LAMMPS:

.. testcode:: python
   :skipif: not HAS_LAMMPS

   from narupatools.ase import ASEDynamics
   from narupatools.lammps import LAMMPSSimulation, atoms_from_lammps_simulation

   # Read in the simulation from an input file and associated data file.
   simulation = LAMMPSSimulation.from_file("in.peptide", "data.peptide")

   # Create an ASE atoms object from information provided by the simulation.
   # This atoms object has a calculator that will use the provided simulation for forces/energies
   atoms = atoms_from_lammps_simulation(simulation)

   # Create a velocity verlet simulation in ASE
   dynamics = ASEDynamics.create_velocity_verlet(atoms, timestep=0.002)

   # Run some dynamics
   dynamics.run(100)

Running Directly
================

As with OpenMM, there is also the option to run a LAMMPS simulation dynamics directly, allowing the use of integrators and other features found in LAMMPS. This is done by creating a :class:`~narupatools.lammps.dynamics.LAMMPSDynamics` object, passing in the simulation. The dynamics object wraps the simulation with all the standard features of the narupatools :class:`~narupatools.imd.dynamics.InteractiveSimulationDynamics`, including throttling the simulation to run at a certain speed, play, pause and applying interactive forces.

Interactive forces are applied to the LAMMPS simulation using the `addforce` fix. These additions are applied automatically without altering the LAMMPS input file.

.. testcode:: python
   :skipif: not HAS_LAMMPS

   from narupatools.ase import ASEDynamics
   from narupatools.lammps import LAMMPSSimulation, LAMMPSDynamics

   # Read in the simulation from an input file and associated data file.
   simulation = LAMMPSSimulation.from_file("in.peptide", "data.peptide")

   # Wrap the simulation in dynamics that can be run.
   dynamics = LAMMPSDynamics(simulation)

   # Run some dynamics
   dynamics.run(100)