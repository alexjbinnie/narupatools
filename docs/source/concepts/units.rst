Units and Conversions
=====================

The units used by *narupatools* are the same as used for Narupa and OpenMM:

.. list-table:: Units of *narupatools*
   :header-rows: 1

   * - Quantity
     - Unit
     - SI Unit
   * - Length
     - Nanometers (:math:`\text{nm}`)
     - :math:`10^{-9} \text{m}`
   * - Time
     - Picoseconds (:math:`\text{ps}`)
     - :math:`10^{-12} \text{s}`
   * - Mass
     - Dalton (:math:`\text{Da}`)
     - :math:`1.661 \times 10^{-27} \text{kg}`
   * - Charge
     - Elementary charge (:math:`\text{e}`)
     - :math:`1.602 \times 10^{-19} \text{C}`
   * - Temperature
     - Kelvin (:math:`\text{K}`)
     -
   * - Angle
     - Radians (:math:`\text{rad}`)
     -
   * - Energy
     - Kilojoules per mole (:math:`\text{kJ} \text{mol}^{-1}`)
     - :math:`6.022 \times 10^{-27} \text{J}`
   * - Force
     - Kilojoules per mole per nanometer (:math:`\text{kJ} \text{mol}^{-1} \text{nm}^{-1}`)
     -
   * - Velocity
     - Nanometers per picosecond (:math:`\text{nm} \text{ps}^{-1}`)
     -
   * - Acceleration
     - Nanometers per picosecond squared (:math:`\text{nm} \text{ps}^{-2}`)
     -

These units are **consistent**, in that the units of mass times acceleration are equal to the units of force. This is not true in certain other packages such as MDAnalysis.

Unit conversions are normally done by either having defined constants (such as :code:`NM_TO_A` to convert nanometers to angstroms), or through wrapping quantities like OpenMM's `Quantity`.

This package provides useful methods for converting between the various unit systems, without having to know exactly what units they use and hence which constants to use. Each package has a :class:`~UnitSystem` that stores the units used for each quantity:

* Narupa/narupatools has the unit system :obj:`~narupatools.core.units.UnitsNarupa`. This uses the unit system as described above.
* OpenMM has the unit system :obj:`~narupatools.openmm.units.UnitsOpenMM`. This is the same unit system as Narupa.
* ASE has the unit system :obj:`~narupatools.ase.units.UnitsASE`. The unit system of ASE is described `here <https://wiki.fysik.dtu.dk/ase/ase/units.html>`_.
* MDAnalysis has the unit system :obj:`~narupatools.ase.units.UnitsMDAnalysis`. It is similar to Narupa/OpenMM, except it uses degrees for angles and angstrom for length.
* MDTraj has the unit system :obj:`~narupatools.ase.units.UnitsMDTraj`. It is similar to Narupa/OpenMM, except it uses degrees for angles.
* LAMMPS actually has several unit systems (as described `here <https://lammps.sandia.gov/doc/units.html>`_. Each of these exists as a separate unit system (such as :obj:`~narupatools.lammps.units.UnitsLAMMPSNano`), but there is also a function :func:`narupatools.lammps.get_unit_system(str)` that converts the LAMMPS name such as 'nano' into the corresponding unit system.

Each of these are instances of a :obj:`~narupatools.core.units.UnitSystem`. A unit system has attributes for all the main quantities that might need to be computed, such as :obj:`~narupatools.core.units.UnitSystem.force` or :obj:`~narupatools.core.units.UnitSystem.length`. The value returned by these will be the value this quantity takes expressed in standard Narupa units. For example, as MDAnalysis uses angstroms, :code:`UnitsMDAnalysis.length` has a value of 0.1.

Conversions between these are done by creating a unit system conversion using the :code:`>>` operator:

.. code-block:: python

    ASEToNarupa = UnitsASE >> UnitsNarupa

    positions_narupa = atoms.get_positions() * ASEToNarupa.length

This object acts in a similar way to the unit system, except now properties such as :code:`.length` would give the conversion factor from lengths in ASE units to lengths in Narupa units.

There are also units defined in :mod:`narupatools.core.units` such as :obj:`~narupatools.core.units.electronvolt`. These units can be multiplied, divided and exponentiated together to get another unit. When multiplied on the left by a float, they represent that value in this units, expressed in standard Narupa units. For example:

.. code-block:: python

    from narupatools.core.units import angstrom

    value_nm = 2.42 * angstrom

Units also override the :code:`>>` operator and hence can be used to get conversion factors in a pythonic way:

.. code-block:: python

    from narupatools.core.units import joule, electronvolt

    joules_to_electronvolts = joule >> electronvolt
