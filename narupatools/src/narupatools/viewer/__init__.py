import importlib

has_lammps = importlib.util.find_spec("pythreejs") is not None

if not has_lammps:
    raise ImportError("narupatools.viewer requires pythreejs to be installed.")

from ._scene import show

__all__ = ["show"]
