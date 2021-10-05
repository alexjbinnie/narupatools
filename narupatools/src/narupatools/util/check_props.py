from typing import List, Union

from narupatools.util.properties import auto


class MyClass:
    @auto
    def my_float(self) -> float:
        """ """

    @auto
    def my_int(self) -> int:
        """ """

    @auto
    def my_str(self) -> str:
        """ """

    @auto
    def my_list(self) -> List[str]:
        """ """

    @auto
    def my_union(self) -> Union[int, str]:
        """ """


obj = MyClass()

obj.my_float = 3.0

obj.my_int = 32

obj.my_str = "abc"

obj.my_list = (2, 4, 7)

obj.my_list = range(5)

print(obj.my_list)
