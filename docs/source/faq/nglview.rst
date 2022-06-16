Using with nglview
==================

The nglview package provides a dynamic widget in a Jupyter notebook for visualising molecular systems.

Narupatools provides handy utilities for using nglview with various constructs.

How do I...
-----------

... visualise a dynamics object?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A dynamics object can be shown using the :obj:`~narpatools.nglview.show_dynamics` function. This generates an nglview
widget that updates its contents as the dynamics object is changed.

.. testsetup:: dynamics

   from ase import Atoms
   from narupatools.ase import ASEDynamics

   dynamics = ASEDynamics.create_velocity_verlet(Atoms(), timestep=0.1)

.. testcode:: dynamics

   from narupatools.nglview import show_dynamics

   show_dynamics(dynamics)

... use shorthand Jupyter magick to use nglview
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The narupatools package provides IPython magick to allow you to quickly render a system. Firstly, this must be setup:

.. code-block:: python

   # Load the narupatools IPython extension
   %load_ext narupatools

Then, instead of having to import the correct ``show_`` function, the ``%ngl`` magick will automatically call the correct
version:

.. code-block:: python

   # Render dynamics using nglview
   %ngl dynamics