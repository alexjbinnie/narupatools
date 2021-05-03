######
OpenMM
######

`OpenMM <http://openmm.org/>`_ is a high-performance simulation engine for molecular dynamics. It is the recommended simulation engine for use with *narupatools*.

As with Narupa, it can be used as either a calculator for ASE or directly.

As an ASE Calculator
====================

The OpenMM simulation can be treated as a calculator for an ASE simulation. An ASE calculator describes how forces and energies are computed. By using OpenMM as a calculator, you are leveraging the multitude of force fields available in OpenMM, whilst using the integrators present in ASE. Whilst slower than using OpenMM directly, as ASE is written in python it is easier to customize and prototype.

The easiest way to construct a runnable molecular dynamics this way is using :meth:`~narupatools.ase.openmm.dynamics.OpenMMASEDynamics.from_simulation`, or :meth:`~narupatools.ase.openmm.dynamics.OpenMMASEDynamics.from_xml_file`:

.. code-block:: python

   from narupatools.ase.openmm import OpenMMASEDynamics

   dynamics = OpenMMASEDynamics.from_xml_file("neuraminidase.xml")

   dynamics.run(100)

By default, this will attempt to use a similar integrator to that defined in the OpenMM xml. If that is not possible, a default langevin integrator will be used.

Running OpenMM directly
=======================

OpenMM can also be run directly, by using :class:`~narupatools.openmm.dynamics.OpenMMDynamics`. This is more performant than running through OpenMM, and allows use of the integrators and other features available in OpenMM.

.. code-block:: python

   from narupatools.openmm import OpenMMDynamics

   dynamics = OpenMMDynamics.from_xml_file("neuraminidase.xml")

   dynamics.run(100)

