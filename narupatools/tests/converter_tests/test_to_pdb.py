from narupa.trajectory import FrameData

from narupatools.frame._converter import frame_to_pdb_string


def test_pdb_empty_frame():
    frame_to_pdb_string(FrameData())


def test_pdb_line_length(frame):
    pdb = frame_to_pdb_string(frame)
    lines = pdb.splitlines()
    for line in lines:
        assert len(line) == 80, f"Line is not of length 80: `{line}`"
