import contextlib
from typing import Type, TypeVar

import numpy as np
import numpy.typing as npt
from rdkit import Chem
from rdkit.Chem import AllChem, CombineMols
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdShapeHelpers import ComputeConfDimsAndOffset
from rdkit.Geometry import Point3D

from narupatools.frame import convert


def _scatter_rectangles(rectangles: npt.ArrayLike) -> np.ndarray:
    """
    Generate positions to arrange a set of rectangles.

    This takes in a set of rectangles of the form [xsize, ysize, zsize]. It then attempts
    to arrange those rectangles such that they are not intersecting.
    """
    rectangles = np.asfarray(rectangles)

    assert len(rectangles.shape) == 2

    count = rectangles.shape[0]
    dim = rectangles.shape[1]

    points = np.zeros((count, dim), dtype=float)

    # Initial bounding box is just the first rectangle.
    bounding_box = np.array([np.zeros(dim), rectangles[0]])

    for i in range(1, count):
        # Size of the cuboid
        size = rectangles[i]

        r = np.random.random_sample(dim) - 0.5
        r /= bounding_box[1]
        r /= np.sqrt(np.max(r ** 2))

        points[i] = (
            bounding_box[0]
            + 0.5 * (bounding_box[1] - size)
            + 0.5 * r * (bounding_box[1] + size)
        )

        bb_corner = np.minimum(bounding_box[0], points[i])
        bounding_box = np.array(
            [
                bb_corner,
                np.maximum(bounding_box[0] + bounding_box[1], points[i] + size)
                - bb_corner,
            ]
        )

    midpoint = bounding_box[0] + 0.5 * bounding_box[1]

    return points - midpoint  # type: ignore


_TType = TypeVar("_TType")


def generate_from_smiles(
    *smiles: str, add_hydrogens: bool = True, output_type: Type[_TType] = Mol
) -> _TType:
    """
    Generate an ASE Atoms object from a SMILES string.

    This leverages RDKit's ability to generate molecules from SMILES
    strings to create a system.

    The structures are minimized using the UFF force field.

    :param smiles: SMILES string describing the system.
    :param add_hydrogens: Should implicit hydrogens be added to the system.
    :return: Atoms object representing system.
    """
    # Generate each MOL separately.
    mols = []
    for index, smile in enumerate(smiles):
        mol = Chem.MolFromSmiles(smile)
        if add_hydrogens:
            mol = AllChem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        # Minimize to avoid clashes
        with contextlib.suppress(RuntimeError):
            AllChem.UFFOptimizeMolecule(
                mol, ignoreInterfragInteractions=True, maxIters=2000
            )

        monomer_info = Chem.AtomPDBResidueInfo()
        monomer_info.SetResidueName("UNL")
        monomer_info.SetResidueNumber(index + 1)
        for atom in mol.GetAtoms():
            atom.SetMonomerInfo(monomer_info)

        mols.append(mol)

    # Space out the molecules through space
    rects = []
    for mol in mols:
        size, corner = ComputeConfDimsAndOffset(mol.GetConformer(0))
        rects.append([corner, size])
    rects = np.array(rects)
    offsets = _scatter_rectangles(rects[:, 1])
    mol = Mol()
    for newmol, offset, rect in zip(mols, offsets, rects):
        translation = offset - rect[0]
        mol = CombineMols(mol, newmol, Point3D(*translation))

    # Minimize the structure of the molecule
    Chem.SanitizeMol(mol)
    with contextlib.suppress(RuntimeError):
        AllChem.UFFOptimizeMolecule(
            mol, ignoreInterfragInteractions=False, maxIters=100
        )

    return convert(mol, output_type)
