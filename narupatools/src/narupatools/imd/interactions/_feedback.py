from narupatools.core.properties import str_property, float_property
from narupatools.state import SharedStateObject


class InteractionFeedback(SharedStateObject):

    @str_property
    def interaction_type(self) -> str:
        """Internal type of the interaction."""

    @float_property
    def total_work(self) -> float:
        """Total work performed by the interaction in kilojoules per mole."""

    @float_property
    def potential_energy(self) -> float:
        """Potential energy of the interaciton in kilojoules per mole."""
