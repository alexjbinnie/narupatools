from narupatools.frame.fields import BondPairs


def test_str_equality():
    assert BondPairs == "bond.pairs"


def test_str_inequality():
    assert BondPairs != "particle.positions"


def test_str_set_membership():
    assert BondPairs in {"bond.pairs"}


def test_key_set_membership():
    assert BondPairs in {BondPairs}


def test_key_str_set_membership():
    assert "bond.pairs" in {BondPairs}


def test_key_dict_key():
    dict = {
        "bond.pairs": 2
    }
    assert dict[BondPairs] == 2