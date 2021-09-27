from narupatools.state import SharedStateObject
from narupatools.util import properties


class InteractionFeedback(SharedStateObject):
    """Feedback for an interaction containing information on its state."""

    @properties.string
    def interaction_type(self) -> str:
        """Internal type of the interaction."""

    @properties.number
    def total_work(self) -> float:
        """Total work performed by the interaction in kilojoules per mole."""

    @properties.number
    def potential_energy(self) -> float:
        """Potential energy of the interaciton in kilojoules per mole."""
