from types import TracebackType
from typing import Any, Dict, List, Optional, Type

from lammps.core import lammps

class System:
    @property
    def units(self) -> str: ...

class OutputCapture:
    def __enter__(self) -> OutputCapture: ...
    def __exit__(
        self,
        t: Optional[Type[BaseException]] = ...,
        value: Optional[BaseException] = ...,
        traceback: Optional[TracebackType] = ...,
    ) -> bool: ...
    @property
    def output(self) -> str: ...

class PyLammps:
    enable_cmd_history: bool

    _cmd_history: List[str]
    @property
    def system(self) -> System: ...
    @property
    def lmp(self) -> lammps: ...
    def command(self, command: str) -> None: ...
    def run(self, steps: int) -> None: ...
    def file(self, filename: str) -> None: ...
    def units(self, units: str) -> None: ...
    def atom_modify(self, keywords: str) -> None: ...
    @property
    def computes(self) -> List[Dict[str, str]]: ...
    def info(self, type: str) -> List[str]: ...
    def eval(self, expression: str) -> Any: ...
    def close(self) -> None: ...