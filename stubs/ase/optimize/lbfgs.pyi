from typing import Optional

from ase import Atoms

class LBFGS:
    def __init__(self, atoms: Atoms, logfile: Optional[str] = ...): ...
    def run(self, fmax: float = ..., steps: Optional[int] = ...): ...
