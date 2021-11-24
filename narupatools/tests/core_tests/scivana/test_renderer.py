from narupatools.app.scivana import Color, Render, Renderer


def test_renderers():
    assert Renderer.create(
        color=Color.by_element(scheme="narupa element"), render=Render.ball_and_stick()
    ).serialize() == {
        "render": {"type": "ball and stick"},
        "color": {"type": "particle element", "scheme": "narupa element"},
    }
