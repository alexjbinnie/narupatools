from __future__ import annotations

from typing import ClassVar, List, Mapping, Optional, Union

import narupatools.util.properties as properties
from narupatools.state import serialize_as
from narupatools.state.typing import Serializable

from ._subgraph import SubgraphObject
from ._typing import Element, Gradient, SingleColor


class Color(SubgraphObject):
    @classmethod
    def by_element(
        cls,
        *,
        scheme: Optional[Union[str, Mapping[Element, SingleColor]]] = None,
        **kwargs: Serializable
    ) -> ColorByElement:
        """
        Color particles by their atomic element.

        :param scheme: Color scheme mapping elements to colors.
        """
        return ColorByElement(scheme=scheme, **kwargs)

    @classmethod
    def by_particle_type(
        cls,
        *,
        scheme: Optional[Union[str, Mapping[str, SingleColor]]] = None,
        **kwargs: Serializable
    ) -> ColorByParticleType:
        """
        Color by the particle type.

        :param scheme: Color scheme mapping particle types to colors.
        """
        return ColorByParticleType(scheme=scheme, **kwargs)

    @classmethod
    def by_residue_name(
        cls,
        *,
        scheme: Optional[Union[str, Mapping[str, SingleColor]]] = None,
        **kwargs: Serializable
    ) -> ColorByResidueName:
        """
        Color by reside name.

        :param scheme: Color scheme mapping residue names to colors.
        """
        return ColorByResidueName(scheme=scheme, **kwargs)

    @classmethod
    def by_secondary_structure(
        cls,
        *,
        scheme: Optional[Union[str, Mapping[str, SingleColor]]] = None,
        **kwargs: Serializable
    ) -> ColorBySecondaryStructure:
        """
        Color by secondary structure assignment.

        :param scheme: Color scheme mapping secondary structure assignments to colors.
        """
        return ColorBySecondaryStructure(scheme=scheme, **kwargs)

    @classmethod
    def by_float(
        cls,
        *,
        input: Optional[List[float]] = None,
        minimum: Optional[float] = None,
        maximum: Optional[float] = None,
        gradient: Optional[Gradient] = None,
        **kwargs: Serializable
    ) -> ColorByFloatGradient:
        """
        Color by an arbitrary field such as particle charge.

        :param input: Key of the property which should provide the data, preceded by '$'.
        :param minimum: Value which should be mapped to the start of the gradient.
        :param maximum: Value which should be mapped to the end of the gradient
        :param gradient: Gradient of colors to use.
        """
        return ColorByFloatGradient(
            input=input, minimum=minimum, maximum=maximum, gradient=gradient, **kwargs
        )

    @classmethod
    def by_particle_index(
        cls, *, gradient: Optional[Gradient] = None, **kwargs: Serializable
    ) -> ColorByParticleIndex:
        """
        Color by particle index.

        :param gradient: Gradient of colors to use.
        """
        return ColorByParticleIndex(gradient=gradient, **kwargs)

    @classmethod
    def by_residue_index(
        cls, *, gradient: Optional[Gradient] = None, **kwargs: Serializable
    ) -> ColorByResidueIndex:
        """
        Color by residue index.

        :param gradient: Gradient of colors to use.
        """
        return ColorByResidueIndex(gradient=gradient, **kwargs)

    @classmethod
    def goodsell(
        cls, *, scheme: Optional[List[SingleColor]] = None, **kwargs: Serializable
    ) -> Goodsell:
        """
        Color in the Goodsell style, colored by chains with carbons being paler.

        :param scheme: List of colors to apply sequentially to entities.
        """
        return Goodsell(scheme=scheme, **kwargs)


class ColorByElement(Color):
    """Color particles by their atomic element."""

    _subgraph_ids: ClassVar[List[str]] = ["particle element", "cpk", "element"]

    def __init__(
        self,
        *,
        scheme: Optional[Union[str, Mapping[Element, SingleColor]]] = None,
        **kwargs: Serializable
    ):
        """
        Color particles by their atomic element.

        :param scheme: Color scheme mapping elements to colors.
        """
        super().__init__(**kwargs)
        if scheme is not None:
            self.scheme = scheme

    @serialize_as("scheme")
    @properties.auto
    def scheme(self) -> Union[str, Mapping[Element, SingleColor]]:
        """Color scheme mapping elements to colors."""


class ColorByParticleType(Color):
    """Color by the particle type."""

    _subgraph_ids: ClassVar[List[str]] = ["particle type", "type"]

    def __init__(
        self,
        *,
        scheme: Optional[Union[str, Mapping[str, SingleColor]]] = None,
        **kwargs: Serializable
    ):
        """
        Color by the particle type.

        :param scheme: Color scheme mapping particle types to colors.
        """
        super().__init__(**kwargs)
        if scheme is not None:
            self.scheme = scheme

    @serialize_as("scheme")
    @properties.auto
    def scheme(self) -> Union[str, Mapping[str, SingleColor]]:
        """Color scheme mapping particle types to colors."""


class ColorByResidueName(Color):
    """Color by reside name."""

    _subgraph_ids: ClassVar[List[str]] = ["residue name", "resname"]

    def __init__(
        self,
        *,
        scheme: Optional[Union[str, Mapping[str, SingleColor]]] = None,
        **kwargs: Serializable
    ):
        """
        Color by reside name.

        :param scheme: Color scheme mapping residue names to colors.
        """
        super().__init__(**kwargs)
        if scheme is not None:
            self.scheme = scheme

    @serialize_as("scheme")
    @properties.auto
    def scheme(self) -> Union[str, Mapping[str, SingleColor]]:
        """Color scheme mapping residue names to colors."""


class ColorBySecondaryStructure(Color):
    """Color by secondary structure assignment."""

    _subgraph_ids: ClassVar[List[str]] = ["secondary structure"]

    def __init__(
        self,
        *,
        scheme: Optional[Union[str, Mapping[str, SingleColor]]] = None,
        **kwargs: Serializable
    ):
        """
        Color by secondary structure assignment.

        :param scheme: Color scheme mapping secondary structure assignments to colors.
        """
        super().__init__(**kwargs)
        if scheme is not None:
            self.scheme = scheme

    @serialize_as("scheme")
    @properties.auto
    def scheme(self) -> Union[str, Mapping[str, SingleColor]]:
        """Color scheme mapping secondary structure assignments to colors."""


class ColorByFloatGradient(Color):
    """Color by an arbitrary field such as particle charge."""

    _subgraph_ids: ClassVar[List[str]] = ["float data gradient"]

    def __init__(
        self,
        *,
        input: Optional[List[float]] = None,
        minimum: Optional[float] = None,
        maximum: Optional[float] = None,
        gradient: Optional[Gradient] = None,
        **kwargs: Serializable
    ):
        """
        Color by an arbitrary field such as particle charge.

        :param input: Key of the property which should provide the data, preceded by '$'.
        :param minimum: Value which should be mapped to the start of the gradient.
        :param maximum: Value which should be mapped to the end of the gradient
        :param gradient: Gradient of colors to use.
        """
        super().__init__(**kwargs)
        if input is not None:
            self.input = input
        if minimum is not None:
            self.minimum = minimum
        if maximum is not None:
            self.maximum = maximum
        if gradient is not None:
            self.gradient = gradient

    @serialize_as("input")
    @properties.auto
    def input(self) -> List[float]:
        """Key of the property which should provide the data, preceded by '$'."""

    @serialize_as("minimum")
    @properties.auto
    def minimum(self) -> float:
        """Value which should be mapped to the start of the gradient."""

    @serialize_as("maximum")
    @properties.auto
    def maximum(self) -> float:
        """Value which should be mapped to the end of the gradient"""

    @serialize_as("gradient")
    @properties.auto
    def gradient(self) -> Gradient:
        """Gradient of colors to use."""


class ColorByParticleIndex(Color):
    """Color by particle index."""

    _subgraph_ids: ClassVar[List[str]] = ["particle index", "particle index in system"]

    def __init__(self, *, gradient: Optional[Gradient] = None, **kwargs: Serializable):
        """
        Color by particle index.

        :param gradient: Gradient of colors to use.
        """
        super().__init__(**kwargs)
        if gradient is not None:
            self.gradient = gradient

    @serialize_as("gradient")
    @properties.auto
    def gradient(self) -> Gradient:
        """Gradient of colors to use."""


class ColorByResidueIndex(Color):
    """Color by residue index."""

    _subgraph_ids: ClassVar[List[str]] = ["residue index", "residue index in system"]

    def __init__(self, *, gradient: Optional[Gradient] = None, **kwargs: Serializable):
        """
        Color by residue index.

        :param gradient: Gradient of colors to use.
        """
        super().__init__(**kwargs)
        if gradient is not None:
            self.gradient = gradient

    @serialize_as("gradient")
    @properties.auto
    def gradient(self) -> Gradient:
        """Gradient of colors to use."""


class Goodsell(Color):
    """Color in the Goodsell style, colored by chains with carbons being paler."""

    _subgraph_ids: ClassVar[List[str]] = ["goodsell"]

    def __init__(
        self, *, scheme: Optional[List[SingleColor]] = None, **kwargs: Serializable
    ):
        """
        Color in the Goodsell style, colored by chains with carbons being paler.

        :param scheme: List of colors to apply sequentially to entities.
        """
        super().__init__(**kwargs)
        if scheme is not None:
            self.scheme = scheme

    @serialize_as("scheme")
    @properties.auto
    def scheme(self) -> List[SingleColor]:
        """List of colors to apply sequentially to entities."""
