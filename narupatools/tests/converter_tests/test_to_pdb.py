from narupatools.frame.converter import frame_to_pdb_string


def test_to_pdb(frame):
    pdb = frame_to_pdb_string(frame)
