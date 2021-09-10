from typing import Optional, Union

import numpy as np
import numpy.typing as npt
from ase import Atoms
from ase.calculators.calculator import PropertyNotImplementedError
from ase.md.md import MolecularDynamics

from narupatools.physics.quaternion import from_vector_part, quaternion
from narupatools.physics.typing import Vector3, Vector3Array
from narupatools.physics.vector import normalized
from ._rotations import set_angular_momenta, get_torques, get_angular_momenta, get_principal_moments, get_rotations, \
    calculate_angular_velocity, set_rotations


def right_multiply(
        v: Union[Vector3, Vector3Array], q: Union[npt.NDArray[quaternion], quaternion], /
) -> npt.NDArray[quaternion]:
    return from_vector_part(v) * q  # type: ignore


class RotationalVelocityVerletIntegrator(MolecularDynamics):
    """
    Symplectic velocity-verlet integrator that supports particle rotations.
    """

    def __init__(self, atoms: Atoms, *, timestep: float):
        super().__init__(atoms=atoms, timestep=timestep, trajectory=None)

        if not self.atoms.has("angmom"):
            set_angular_momenta(atoms, np.zeros([len(self.atoms), 3]))

    def step(
            self,
            forces: Optional[Vector3Array] = None,
            torques: Optional[Vector3Array] = None,
    ) -> None:

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

        omega = calculate_angular_velocity(
            principal_moments=I, angular_momenta=L, orientations=q
        )

        dqdt = 0.5 * right_multiply(omega, q)  # type: ignore [operator]

        qfull = normalized(q + self.dt * dqdt)
        qhalf = normalized(q + 0.5 * self.dt * dqdt)

        omega = calculate_angular_velocity(
            principal_moments=I, angular_momenta=L, orientations=qhalf
        )

        dqdt = 0.5 * right_multiply(omega, q)  # type: ignore [operator]

        qhalf = normalized(qhalf + 0.5 * self.dt * dqdt)

        q = normalized(2 * qhalf - qfull)

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
