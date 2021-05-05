import pytest
from infinite_sets import InfiniteSet, everything
from narupa.trajectory import FrameData

from narupatools.core.trajectory import TrajectoryPlayback

ALL_FIELDS = everything()


class TestTrajectory(TrajectoryPlayback):
    def __init__(self):
        super().__init__(looping=True)

    def _trajectory_length(self) -> int:
        return 5

    def get_frame(self, fields: InfiniteSet[str] = ALL_FIELDS) -> FrameData:
        data = FrameData()
        data.values["index"] = self.index
        return data


@pytest.fixture
def trajectory():
    return TestTrajectory()


def test_initial_index(trajectory):
    assert trajectory.index == 0
    assert trajectory.get_frame().values["index"] == 0


def test_set_index(trajectory):
    trajectory.index = 3
    assert trajectory.index == 3
    assert trajectory.get_frame().values["index"] == 3


def test_set_index_out_of_range(trajectory):
    with pytest.raises(IndexError):
        trajectory.index = 14


def test_set_index_negative(trajectory):
    with pytest.raises(IndexError):
        trajectory.index = -11


def test_increment_step(trajectory):
    trajectory.step()
    assert trajectory.index == 1
    assert trajectory.get_frame().values["index"] == 1


def test_increment_step_multiple(trajectory):
    for _ in range(4):
        trajectory.step()
    assert trajectory.index == 4
    assert trajectory.get_frame().values["index"] == 4


def test_looping(trajectory):
    for _ in range(7):
        trajectory.step()
    assert trajectory.index == 2
    assert trajectory.get_frame().values["index"] == 2


def test_non_looping(trajectory):
    trajectory.looping = False
    for _ in range(7):
        trajectory.step()
    assert trajectory.index == 4
    assert trajectory.get_frame().values["index"] == 4


def test_restart(trajectory):
    for _ in range(4):
        trajectory.step()
    trajectory.restart()
    assert trajectory.index == 0
    assert trajectory.get_frame().values["index"] == 0
