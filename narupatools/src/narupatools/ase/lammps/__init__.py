from ase import Atoms

from ._calculator import LAMMPSCalculator
from narupatools.frame import convert
from narupatools.lammps import LAMMPSSimulation


def atoms_from_lammps_simulation(simulation: LAMMPSSimulation) -> Atoms:
    """
    Create an ASE Atoms object based on a LAMMPS simulation.

    :param simulation: LAMMPS simulation to use.
    :return: ASE atoms object with a calculator that references the given simulation.
    """
    atoms = convert(simulation, Atoms)
    calc = LAMMPSCalculator(simulation, atoms)
    atoms.set_calculator(calc)
    return atoms
