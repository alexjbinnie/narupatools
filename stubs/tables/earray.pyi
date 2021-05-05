from typing import Any

from tables.node import Node

class EArray(Node):
    def append(self, sequence: Any) -> None: ...
