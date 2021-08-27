from typing import Collection, Optional

from ase import Atoms

from .calculator import Calculator

class AverageCalculator(Calculator):
    def __init__(self, calcs: Collection[Calculator], atoms: Optional[Atoms] = ...): ...
