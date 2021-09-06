import numpy as np
from ase import Atoms
from ase.calculators.calculator import PropertyNotImplementedError
from ase.md.md import MolecularDynamics
from narupatools.physics.quaternion import quaternion, as_quat_array, from_vector_part

from narupatools.physics.typing import Vector3Array

ANGMOM_ARRAY = "angmom"
ORIENTATION_ARRAY = "rotations"
PRINCIPAL_MOMENTS_ARRAY = "principal_moments"
TORQUES_PROPERTY = "torques"


def set_angular_momenta(atoms: Atoms, momenta: Vector3Array, apply_constraint=True):
    if apply_constraint and len(atoms.constraints) > 0 and momenta is not None:
        momenta = np.array(momenta)
        for constraint in atoms.constraints:
            if hasattr(constraint, 'adjust_angular_momenta'):
                constraint.adjust_angular_momenta(atoms, momenta)
    atoms.set_array(ANGMOM_ARRAY, momenta, float, (3,))


def get_principal_moments(atoms: Atoms):
    if PRINCIPAL_MOMENTS_ARRAY in atoms.arrays:
        return atoms.arrays[PRINCIPAL_MOMENTS_ARRAY].copy()
    else:
        return np.zeros(len(atoms))


def set_principal_moments(atoms: Atoms, value):
    atoms.arrays[PRINCIPAL_MOMENTS_ARRAY] = value


def get_angular_momenta(atoms: Atoms):
    if ANGMOM_ARRAY in atoms.arrays:
        return atoms.arrays[ANGMOM_ARRAY].copy()
    else:
        return np.zeros((len(atoms), 3))


def set_rotations(atoms: Atoms, q: np.ndarray):
    atoms.arrays[ORIENTATION_ARRAY] = q


def get_rotations(atoms: Atoms):
    if ORIENTATION_ARRAY in atoms.arrays:
        array = atoms.arrays[ORIENTATION_ARRAY].copy()
        if array.dtype == quaternion:
            return array
        return as_quat_array(array)
    else:
        return np.repeat(quaternion(1, 0, 0, 0), len(atoms))


def get_torques(atoms: Atoms, apply_constraint=True, md=False):
    if atoms.calc is None:
        raise RuntimeError('Atoms object has no calculator.')

    try:
        torques = atoms.calc.get_property(TORQUES_PROPERTY, atoms)
    except PropertyNotImplementedError:
        torques = np.zeros((len(atoms), 3))
    if apply_constraint:
        for constraint in atoms.constraints:
            if not md or hasattr(constraint, 'adjust_torques'):
                constraint.adjust_torques(atoms, torques)
    return torques


def calculate_angular_velocity(*, principal_moments, angular_momentum, orientation):
    if len(principal_moments.shape) == 1:
        return np.nan_to_num(angular_momentum / principal_moments)  # CHECK!
    else:
        # todo
        raise ValueError

def set_angular_velocities(atoms: Atoms,  angular_velocities):
    principal_moments = get_principal_moments(atoms)
    if len(principal_moments.shape) == 1:
        set_angular_momenta(atoms, angular_velocities * principal_moments)
    else:
        raise ValueError



def right_multiply(v, q):
    return from_vector_part(v) * q


class RotationalVelocityVerletIntegrator(MolecularDynamics):
    """
    Symplectic velocity-verlet integrator that supports particle rotations.
    """

    def __init__(self, atoms: Atoms, *, timestep: float):
        super().__init__(atoms=atoms, timestep=timestep, trajectory=None)

        if not self.atoms.has('angmom'):
            set_angular_momenta(atoms, np.zeros([len(self.atoms), 3]))

    def step(self, forces=None, torques=None):

        atoms = self.atoms

        if forces is None:
            forces = atoms.get_forces(md=True)

        if torques is None:
            torques = get_torques(atoms, md=True)

        p = atoms.get_momenta()
        L = get_angular_momenta(atoms)

        p += 0.5 * self.dt * forces
        L += 0.5 * self.dt * torques

        q = get_rotations(atoms)
        I = get_principal_moments(atoms)

        omega = calculate_angular_velocity(principal_moments=I, angular_momentum=L, orientation=q)

        dqdt = 0.5 * right_multiply(omega, q)

        qfull = np.normalized(q + self.dt * dqdt)
        qhalf = np.normalized(q + 0.5 * self.dt * dqdt)

        omega = calculate_angular_velocity(principal_moments=I, angular_momentum=L, orientation=qhalf)

        dqdt = 0.5 * right_multiply(omega, q)

        qhalf = np.normalized(qhalf + 0.5 * self.dt * dqdt)

        q = np.normalized(2 * qhalf - qfull)

        m = atoms.get_masses()[:, np.newaxis]
        r = atoms.get_positions()

        # if we have constraints then this will do the first part of the
        # RATTLE algorithm:
        atoms.set_positions(r + self.dt * p / m)
        if atoms.constraints:
            p = (atoms.get_positions() - r) * m / self.dt

        set_rotations(atoms, q)

        # Set the momenta in case the force depends on it
        atoms.set_momenta(p, apply_constraint=False)
        set_angular_momenta(atoms, L, apply_constraint=False)

        forces = atoms.get_forces(md=True)
        try:
            torques = get_torques(atoms, md=True)
        except PropertyNotImplementedError:
            torques = np.zeros((len(atoms), 3))

        # Second part of RATTLE will be done here:
        p = atoms.get_momenta()
        L = get_angular_momenta(atoms)

        p += 0.5 * self.dt * forces
        L += 0.5 * self.dt * torques

        atoms.set_momenta(p)
        set_angular_momenta(atoms, L)

        return forces
