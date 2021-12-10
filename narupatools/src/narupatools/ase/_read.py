import ase.io

from narupatools.frame import Reader


class ASEReader(Reader):

    def load_frame(self, *filenames: str, index: int = 0):
        if len(filenames) == 1:
            return ase.io.read(filenames[0], index=index)
        return super().load_frame(*filenames)