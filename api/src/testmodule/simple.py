from typing import overload, Union
import sys

from math import *

import numpy as np

global_nodoc_noannot = 5

global_documented: int = 3
"""My Source"""

if False:
    global_in_if = 9

setattr(sys.modules[__name__], "global_setattr", 2)

globals()["global_setitem_globals"] = 3

class SimpleClass:

    class_var: int = 3
    """My Docstring"""

    instance_annot: str

    def __init__(self):
        self.instance_var = None
        self.instance_var_annot: int = 3
        self.instance_var_doc = 8
        """I have a doc"""

    @property
    def readonly_prop(self):
        pass

    @property
    def prop(self):
        pass

    @prop.setter
    def prop(self, value):
        pass

    @overload
    def overloaded(self, a: int) -> int:
        pass

    @overload
    def overloaded(self, a: int) -> int:
        pass

    def overloaded(self, a: Union[int, float]) -> Union[int, float]:
        ...
