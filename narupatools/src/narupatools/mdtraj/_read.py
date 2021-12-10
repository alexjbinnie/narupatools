from typing import Optional

import numpy as np
from infinite_sets import InfiniteSet
from mdtraj import load_frame, load
from narupa.trajectory import FrameData

from narupatools.frame import convert, TrajectorySource, Reader


class MDTrajReader(Reader):

    def load_frame(self, *filenames: str, index: int = 0, subset: Optional[np.ndarray] = None, fields: InfiniteSet[str] = None):
        if len(filenames) == 1:
            filename = filenames[0]
            topology: Optional[str] = None
        elif len(filenames) == 2:
            filename = filenames[0]
            topology = filenames[1]
        else:
            return super().load_frame(*filenames)
        traj = load_frame(filename, index, top=topology, atom_indices=subset)
        return convert(traj, FrameData, fields=fields)

    def load_traj(self, *filenames: str, subset: Optional[np.ndarray] = None):
        traj = load(list(*filenames), atom_indices=subset)
        return TrajectorySource.create_from_object(traj)