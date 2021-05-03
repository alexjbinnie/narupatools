from typing import Union

import numpy as np

class numpy_wrapper:
    def extract_compute(
        self, cid: str, style: int, type: int
    ) -> Union[float, np.ndarray, None]: ...
