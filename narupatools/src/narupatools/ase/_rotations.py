import numpy as np
import numpy.typing as npt
from ase.atoms import Atoms
from ase.calculators.calculator import PropertyNotImplementedError

from narupatools.physics import as_quaternion_array, quaternion
from narupatools.physics.typing import Vector3Array, Vector3Like

ANGMOM_ARRAY = "angmom"
ORIENTATION_ARRAY = "rotations"
PRINCIPAL_MOMENTS_ARRAY = "principal_moments"
TORQUES_PROPERTY = "torques"


def get_angular_momenta(atoms: Atoms, /) -> Vector3Array:
    """Get the angular momenta from an ASE atoms object in ASE units."""
    try:
        return atoms.get_array(ANGMOM_ARRAY)
    except KeyError:
        return np.zeros((len(atoms), 3))


Atoms.get_angular_momenta = get_angular_momenta  # type: ignore[attr-defined]


def set_angular_momenta(
    atoms: Atoms, momenta: Vector3Array, /, apply_constraint: bool = True
) -> None:
    """
    Set the angular momenta of an ASE atoms object.

    :param atoms: Atoms object to change.
    :param momenta: Angular momenta in ASE units.
    :param apply_constraint: Should constraints by applied.
    """
    if apply_constraint and len(atoms.constraints) > 0 and momenta is not None:
        momenta = np.array(momenta)
        for constraint in atoms.constraints:
            if hasattr(constraint, "adjust_angular_momenta"):
                constraint.adjust_angular_momenta(atoms, momenta)
    atoms.set_array(ANGMOM_ARRAY, momenta, float, (3,))


Atoms.set_angular_momenta = set_angular_momenta  # type: ignore[attr-defined]


def get_principal_moments(atoms: Atoms, /) -> Vector3Array:
    """Get the principal moments of inertia for an ASE atoms object."""
    try:
        return atoms.get_array(PRINCIPAL_MOMENTS_ARRAY)
    except KeyError:
        return np.zeros(len(atoms))


Atoms.get_principal_moments = get_principal_moments  # type: ignore[attr-defined]


def set_principal_moments(atoms: Atoms, value: Vector3Array, /) -> None:
    """
    Set the principal moments of inertia for an ASE atoms object.

    :param atoms: Atoms object to change.
    :param value: New principal moments to set.
    """
    atoms.set_array(PRINCIPAL_MOMENTS_ARRAY, value, dtype=float)


Atoms.set_principal_moments = set_principal_moments  # type: ignore[attr-defined]


def get_rotations(atoms: Atoms, /) -> npt.NDArray[quaternion]:
    """Get the rotations of an ASE atoms object."""
    if ORIENTATION_ARRAY in atoms.arrays:
        array = atoms.arrays[ORIENTATION_ARRAY].copy()
        if array.dtype == quaternion:
            return array  # type: ignore
        return as_quaternion_array(array)
    else:
        return np.repeat(quaternion(1, 0, 0, 0), len(atoms))


Atoms.get_rotations = get_rotations  # type: ignore[attr-defined]


def set_rotations(atoms: Atoms, q: npt.NDArray[quaternion], /) -> None:
    """Set the orientations of each atoms object."""
    atoms.set_array(ORIENTATION_ARRAY, q, dtype=quaternion)


Atoms.set_rotations = set_rotations  # type: ignore[attr-defined]


def get_torques(
    atoms: Atoms, apply_constraint: bool = True, md: bool = False
) -> npt.NDArray[np.float_]:
    """Get the torques for an ASE atoms object."""
    if atoms.calc is None:
        raise RuntimeError("Atoms object has no calculator.")

    try:
        torques = atoms.calc.get_property(TORQUES_PROPERTY, atoms)
    except PropertyNotImplementedError:
        torques = np.zeros((len(atoms), 3))
    if apply_constraint:
        for constraint in atoms.constraints:
            if not md or hasattr(constraint, "adjust_torques"):
                constraint.adjust_torques(atoms, torques)
    return torques  # type: ignore


Atoms.get_torques = get_torques  # type: ignore[attr-defined]


def get_angular_velocities(atoms: Atoms, /) -> Vector3Array:
    """
    Get angular velocities from an ASE atoms object.

    This is the per-particle angular velocites about their individual center of masses.
    """
    return calculate_angular_velocity(
        principal_moments=get_principal_moments(atoms),
        angular_momenta=get_angular_momenta(atoms),
        orientations=get_rotations(atoms),
    )


Atoms.get_angular_velocities = get_angular_velocities  # type: ignore[attr-defined]


def calculate_angular_velocity(
    *,
    principal_moments: Vector3Array,
    angular_momenta: Vector3Array,
    orientations: npt.NDArray[quaternion],
) -> Vector3Array:
    """Calculate per-particle angular velocity."""
    if len(principal_moments.shape) == 1:
        return np.nan_to_num(angular_momenta / principal_moments)  # type: ignore
    else:
        # todo
        raise ValueError


def set_angular_velocities(atoms: Atoms, angular_velocities: Vector3Like, /) -> None:
    """Set per-particle angular velocities."""
    principal_moments = get_principal_moments(atoms)
    if len(principal_moments.shape) == 1:
        set_angular_momenta(atoms, angular_velocities * principal_moments)
    else:
        raise ValueError


Atoms.set_angular_velocities = set_angular_velocities  # type: ignore[attr-defined]
