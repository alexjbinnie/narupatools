# This file is part of narupatools (https://github.com/alexjbinnie/narupatools).
# Copyright (c) Alex Jamieson-Binnie. All rights reserved.
#
# narupatools is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# narupatools is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with narupatools.  If not, see <http://www.gnu.org/licenses/>.

from typing import Any, Collection, Dict, List, Optional

from ase.atoms import Atoms

all_changes: List[str]

class CalculatorSetupError(RuntimeError): ...
class PropertyNotImplementedError(NotImplementedError): ...

class Parameters(dict):
    def __getattr__(self, name: str) -> Any: ...

class Calculator:
    results: Dict[str, Any]
    atoms: Optional[Atoms]
    parameters: Parameters
    def __init__(self, **kwargs: Any): ...
    def calculate(
        self,
        atoms: Optional[Atoms] = ...,
        properties: Collection[str] = ...,
        system_changes: List[str] = ...,
    ) -> None: ...
    def reset(self) -> None: ...
    def get_property(
        self, name: str, atoms: Optional[Atoms] = ..., allow_calculation: bool = ...
    ): ...
