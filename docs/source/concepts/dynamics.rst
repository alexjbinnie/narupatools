Dynamics
========

Molecular dynamics as seen by *narupatools* is a method where a system of particles can be simulated by applying forces and integrating through time. This is a brief overview of some of the concepts involved.

There are different representations of dynamics in *narupatools*, such as :class:`~narupatools.ase.ASEDynamics`. All of these are children of the base :class:`~narupatools.core.dynamics.SimulationDynamics`.

When the dynamics is run, either in a background thread or for a certain number of steps on the main thread, the following occurs:

.. graphviz:: dynamics.dot