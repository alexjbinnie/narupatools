Using with OpenMM
=================

RDKit integration in narupatools provides a couple of ASE forcefields, as well as access to the powerful tools allowing generating systems.

How do I...
-----------

... generate a system from one or more SMILES strings?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Systems can be generated from one or more SMILES strings using the :obj:`narupatools.rdkit.generate_from_smiles` functionality:

.. testcode:: from_smiles

   from narupatools.rdkit import generate_from_smiles

   # Create an RDKit Mol object from the SMILES string for ethane
   mol = generate_from_smiles("CC")

.. testcode:: from_smiles

   from narupa.trajectory import FrameData

   # Create a Narupa FrameData from the SMILES string for two methanols
   frame = generate_from_smiles("CO", "CO", output_type=FrameData)

.. testcode:: from_smiles

   from MDAnalysis import Universe

   # Create an MDAnalysis universe from the SMILES string for ethylene without hydrogens
   universe = generate_from_smiles("C=C", output_type=Universe, add_hydrogens=False)

... run ASE dynamics with the MMFF94 forcefield?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :obj:`narupatools.ase.rdkit.MMFF94Calculator` allows the MMFF94 forcefield as implemented by RDKit to be used to provide forces and energies:

.. testsetup:: mmff94

   from ase import Atoms

   atoms = Atoms()

.. testcode:: mmff94

   from narupatools.ase.rdkit import MMFF94Calculator

   atoms.calc = MMFF94Calculator()

... run ASE dynamics with the UFF forcefield?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :obj:`narupatools.ase.rdkit.UFFCalculator` allows the UFF forcefield to be used with ASE:

.. testsetup:: uff

   from ase import Atoms

   atoms = Atoms()

.. testcode:: uff

   from narupatools.ase.rdkit import UFFCalculator

   atoms.calc = UFFCalculator()

... ignore intermolecular forces with MMFF94 or UFF?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Both the MMFF94 and UFF calculators can ignore intermolecular forces by setting *include_interatomic* to False:

.. testsetup:: nonbonded

   from ase import Atoms

   atoms = Atoms()

.. testcode:: nonbonded

   from narupatools.ase.rdkit import UFFCalculator

   atoms.calc = UFFCalculator(include_interatomic=False)

... set a cutoff for nonbonded forces with MMFF94 or UFF?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Both forcefields support having a cutoff specified (in nanometers):

.. testsetup:: cutoff

   from ase import Atoms

   atoms = Atoms()

.. testcode:: cutoff

   from narupatools.ase.rdkit import MMFF94Calculator

   atoms.calc = MMFF94Calculator(nonbonded_cutoff=2.0)