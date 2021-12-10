from MDAnalysis import Universe
from narupa.trajectory import FrameData

from narupatools.frame import convert, TrajectorySource
from narupatools.frame import Reader


class MDAnalysisReader(Reader):

    def load_frame(self, *filenames: str):
        u = Universe(*filenames)
        return convert(u, FrameData)

    def load_trajectory(self, *filenames: str):
        u = Universe(*filenames)
        return TrajectorySource.create_from_object(u)