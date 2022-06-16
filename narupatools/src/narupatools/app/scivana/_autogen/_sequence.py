# flake8: noqa

from __future__ import annotations

from typing import ClassVar, List

from narupatools.state.typing import Serializable

from ._subgraph import SubgraphObject


class Sequence(SubgraphObject):
    @classmethod
    def protein_alpha_carbons(
        cls, **kwargs: Serializable
    ) -> SequenceByPolypeptideAlphaCarbons:
        """
        Generate sequences through contiguous alpha carbons of polypeptides.

        """
        return SequenceByPolypeptideAlphaCarbons(**kwargs)


class SequenceByPolypeptideAlphaCarbons(Sequence):
    """Generate sequences through contiguous alpha carbons of polypeptides."""

    _subgraph_ids: ClassVar[List[str]] = ["polypeptide", "protein"]

    def __init__(self, **kwargs: Serializable):
        """
        Generate sequences through contiguous alpha carbons of polypeptides.

        """
        super().__init__(**kwargs)
