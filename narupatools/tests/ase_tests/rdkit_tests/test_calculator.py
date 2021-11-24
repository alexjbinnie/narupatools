import pytest
from ase import Atoms

from narupatools.ase import ASEDynamics

pytest.importorskip("rdkit")

from narupatools.ase.rdkit import MMFF94Calculator, UFFCalculator
from narupatools.rdkit import generate_from_smiles


def call_count(f):
    def wrapped(*args, **kwargs):
        wrapped.call_count += 1
        return f(*args, **kwargs)

    wrapped.call_count = 0
    return wrapped


@pytest.fixture(params=[MMFF94Calculator, UFFCalculator])
def calculator(request):
    return request.param


def test_forces(calculator):
    atoms = generate_from_smiles("C", output_type=Atoms)
    atoms.calc = calculator()
    assert atoms.get_forces().shape == (5, 3)


def test_calculation_count(calculator):
    atoms = generate_from_smiles("C", output_type=Atoms)
    atoms.calc = calculator()
    atoms.calc._calculate_forcefield = call_count(atoms.calc._calculate_forcefield)
    _ = atoms.get_potential_energy()
    assert atoms.calc._calculate_forcefield.call_count == 1
    _ = atoms.get_potential_energy()
    assert atoms.calc._calculate_forcefield.call_count == 1
    atoms.positions[3] += [0.0, 0.1, 0.0]
    _ = atoms.get_potential_energy()
    assert atoms.calc._calculate_forcefield.call_count == 2


def test_dynamics_call_count(calculator):
    atoms = generate_from_smiles("C", output_type=Atoms)
    atoms.calc = calculator()
    atoms.calc._calculate_forcefield = call_count(atoms.calc._calculate_forcefield)
    dynamics = ASEDynamics.create_velocity_verlet(atoms, timestep=0.001)
    dynamics.run(100)
    assert atoms.calc._calculate_forcefield.call_count == 101
