#####
RDKit
#####

If `RDKit <https://www.rdkit.org/>`_ is installed, then it can be leveraged as both a forcefield and as a tool to generate a system. This makes it ideal for generating quick systems to test.

.. testcode::

   from narupatools.ase import ASEDynamics
   from narupatools.ase.rdkit import atoms_from_smiles, MMFF94Calculator

   # Generate an ASE atoms object of benzene, based on a SMILES string.
   # Unless specified, hydrogens are added automatically
   benzene = atoms_from_smiles("c1ccccc1")

   # Both MMFF94 and UFF are available as calculators
   benzene.calc = MMFF94Calculator()

   dynamics = ASEDynamics.create_langevin(benzene, timestep=0.1, friction=0.1, temperature=273)