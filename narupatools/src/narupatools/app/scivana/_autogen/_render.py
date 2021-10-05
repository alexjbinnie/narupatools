from __future__ import annotations

from typing import ClassVar, List, Optional

import narupatools.util.properties as properties
from narupatools.state import serialize_as
from narupatools.state.typing import Serializable

from ._subgraph import SubgraphObject
from ._typing import SingleColor


class Render(SubgraphObject):
    @classmethod
    def ball_and_stick(
        cls,
        *,
        particle_scale: Optional[float] = None,
        bond_scale: Optional[float] = None,
        scale: Optional[float] = None,
        color: Optional[SingleColor] = None,
        opacity: Optional[float] = None,
        edge_sharpness: Optional[float] = None,
        **kwargs: Serializable
    ) -> RenderBallAndStick:
        """
        Classic ball and stick representation.

        :param particle_scale: Scaling factor to apply to the balls.
        :param bond_scale: Scaling factor to apply to the bonds.
        :param scale:
        :param color:
        :param opacity:
        :param edge_sharpness:
        """
        return RenderBallAndStick(
            particle_scale=particle_scale,
            bond_scale=bond_scale,
            scale=scale,
            color=color,
            opacity=opacity,
            edge_sharpness=edge_sharpness,
            **kwargs
        )

    @classmethod
    def cycles(
        cls,
        *,
        scale: Optional[float] = None,
        cycle_scale: Optional[float] = None,
        color: Optional[SingleColor] = None,
        opacity: Optional[float] = None,
        bond_scale: Optional[float] = None,
        edge_sharpness: Optional[float] = None,
        **kwargs: Serializable
    ) -> RenderCycles:
        """
        Liquorice renderer with solid filled-in cycles.

        :param scale:
        :param cycle_scale:
        :param color:
        :param opacity:
        :param bond_scale:
        :param edge_sharpness:
        """
        return RenderCycles(
            scale=scale,
            cycle_scale=cycle_scale,
            color=color,
            opacity=opacity,
            bond_scale=bond_scale,
            edge_sharpness=edge_sharpness,
            **kwargs
        )

    @classmethod
    def ribbon(
        cls,
        *,
        particle_widths: Optional[List[float]] = None,
        scale: Optional[float] = None,
        width: Optional[float] = None,
        speed: Optional[float] = None,
        input: Optional[List[SingleColor]] = None,
        color: Optional[SingleColor] = None,
        opacity: Optional[float] = None,
        edge_sharpness: Optional[float] = None,
        **kwargs: Serializable
    ) -> RenderRibbon:
        """


        :param particle_widths:
        :param scale:
        :param width:
        :param speed:
        :param input:
        :param color:
        :param opacity:
        :param edge_sharpness:
        """
        return RenderRibbon(
            particle_widths=particle_widths,
            scale=scale,
            width=width,
            speed=speed,
            input=input,
            color=color,
            opacity=opacity,
            edge_sharpness=edge_sharpness,
            **kwargs
        )

    @classmethod
    def geometric_ribbon(
        cls,
        *,
        particle_widths: Optional[List[float]] = None,
        scale: Optional[float] = None,
        width: Optional[float] = None,
        color: Optional[SingleColor] = None,
        opacity: Optional[float] = None,
        bond_scale: Optional[float] = None,
        edge_sharpness: Optional[float] = None,
        **kwargs: Serializable
    ) -> RenderGeometricRibbon:
        """


        :param particle_widths:
        :param scale:
        :param width:
        :param color:
        :param opacity:
        :param bond_scale:
        :param edge_sharpness:
        """
        return RenderGeometricRibbon(
            particle_widths=particle_widths,
            scale=scale,
            width=width,
            color=color,
            opacity=opacity,
            bond_scale=bond_scale,
            edge_sharpness=edge_sharpness,
            **kwargs
        )

    @classmethod
    def goodsell(
        cls,
        *,
        scale: Optional[float] = None,
        color: Optional[SingleColor] = None,
        opacity: Optional[float] = None,
        **kwargs: Serializable
    ) -> RenderGoodsell:
        """
        Goodsell style renderer which draws unshaded spheres with outlines between distant residues.

        :param scale:
        :param color:
        :param opacity:
        """
        return RenderGoodsell(scale=scale, color=color, opacity=opacity, **kwargs)

    @classmethod
    def hydrogen_caps(
        cls,
        *,
        acceptor_scale: Optional[float] = None,
        donor_scale: Optional[float] = None,
        acceptor_color: Optional[SingleColor] = None,
        donor_color: Optional[SingleColor] = None,
        scale: Optional[float] = None,
        color: Optional[SingleColor] = None,
        opacity: Optional[float] = None,
        **kwargs: Serializable
    ) -> RenderHydrogenBondCaps:
        """
        Render cones to indicate directions which hydrogen bonds should exist.

        :param acceptor_scale:
        :param donor_scale:
        :param acceptor_color:
        :param donor_color:
        :param scale:
        :param color:
        :param opacity:
        """
        return RenderHydrogenBondCaps(
            acceptor_scale=acceptor_scale,
            donor_scale=donor_scale,
            acceptor_color=acceptor_color,
            donor_color=donor_color,
            scale=scale,
            color=color,
            opacity=opacity,
            **kwargs
        )

    @classmethod
    def hyperballs(
        cls,
        *,
        scale: Optional[float] = None,
        tension: Optional[float] = None,
        color: Optional[SingleColor] = None,
        opacity: Optional[float] = None,
        bond_scale: Optional[float] = None,
        edge_sharpness: Optional[float] = None,
        **kwargs: Serializable
    ) -> RenderHyperballs:
        """
        Draws spheres connected by hyperboloids.

        :param scale:
        :param tension:
        :param color:
        :param opacity:
        :param bond_scale:
        :param edge_sharpness:
        """
        return RenderHyperballs(
            scale=scale,
            tension=tension,
            color=color,
            opacity=opacity,
            bond_scale=bond_scale,
            edge_sharpness=edge_sharpness,
            **kwargs
        )

    @classmethod
    def liquorice(
        cls,
        *,
        scale: Optional[float] = None,
        color: Optional[SingleColor] = None,
        opacity: Optional[float] = None,
        bond_scale: Optional[float] = None,
        edge_sharpness: Optional[float] = None,
        **kwargs: Serializable
    ) -> RenderLiquorice:
        """
        Render stick representation with rounded bonds.

        :param scale:
        :param color:
        :param opacity:
        :param bond_scale:
        :param edge_sharpness:
        """
        return RenderLiquorice(
            scale=scale,
            color=color,
            opacity=opacity,
            bond_scale=bond_scale,
            edge_sharpness=edge_sharpness,
            **kwargs
        )

    @classmethod
    def noodles(
        cls,
        *,
        scale: Optional[float] = None,
        color: Optional[SingleColor] = None,
        opacity: Optional[float] = None,
        edge_sharpness: Optional[float] = None,
        **kwargs: Serializable
    ) -> RenderNoodles:
        """
        Liquorice-like renderer with bonds curved based on neighbouring bonds.

        :param scale:
        :param color:
        :param opacity:
        :param edge_sharpness:
        """
        return RenderNoodles(
            scale=scale,
            color=color,
            opacity=opacity,
            edge_sharpness=edge_sharpness,
            **kwargs
        )

    @classmethod
    def peptide_planes(
        cls,
        *,
        scale: Optional[float] = None,
        color: Optional[SingleColor] = None,
        opacity: Optional[float] = None,
        bond_scale: Optional[float] = None,
        edge_sharpness: Optional[float] = None,
        **kwargs: Serializable
    ) -> RenderPeptidePlanes:
        """
        Render solid planes aligned with the peptide bonds.

        :param scale:
        :param color:
        :param opacity:
        :param bond_scale:
        :param edge_sharpness:
        """
        return RenderPeptidePlanes(
            scale=scale,
            color=color,
            opacity=opacity,
            bond_scale=bond_scale,
            edge_sharpness=edge_sharpness,
            **kwargs
        )

    @classmethod
    def spheres(
        cls,
        *,
        particle_scale: Optional[float] = None,
        outline_depth: Optional[float] = None,
        outline_width: Optional[float] = None,
        scale: Optional[float] = None,
        color: Optional[SingleColor] = None,
        opacity: Optional[float] = None,
        **kwargs: Serializable
    ) -> RenderSpheres:
        """
        Render spheres without bonds.

        :param particle_scale:
        :param outline_depth:
        :param outline_width:
        :param scale:
        :param color:
        :param opacity:
        """
        return RenderSpheres(
            particle_scale=particle_scale,
            outline_depth=outline_depth,
            outline_width=outline_width,
            scale=scale,
            color=color,
            opacity=opacity,
            **kwargs
        )

    @classmethod
    def tube(
        cls,
        *,
        scale: Optional[float] = None,
        color: Optional[SingleColor] = None,
        opacity: Optional[float] = None,
        edge_sharpness: Optional[float] = None,
        **kwargs: Serializable
    ) -> RenderTube:
        """
        Render solid curved tube through alpha carbons.

        :param scale:
        :param color:
        :param opacity:
        :param edge_sharpness:
        """
        return RenderTube(
            scale=scale,
            color=color,
            opacity=opacity,
            edge_sharpness=edge_sharpness,
            **kwargs
        )


class RenderBallAndStick(Render):
    """Classic ball and stick representation."""

    _subgraph_ids: ClassVar[List[str]] = ["ball and stick", "ball & stick"]

    def __init__(
        self,
        *,
        particle_scale: Optional[float] = None,
        bond_scale: Optional[float] = None,
        scale: Optional[float] = None,
        color: Optional[SingleColor] = None,
        opacity: Optional[float] = None,
        edge_sharpness: Optional[float] = None,
        **kwargs: Serializable
    ):
        """
        Classic ball and stick representation.

        :param particle_scale: Scaling factor to apply to the balls.
        :param bond_scale: Scaling factor to apply to the bonds.
        :param scale:
        :param color:
        :param opacity:
        :param edge_sharpness:
        """
        super().__init__(**kwargs)
        if particle_scale is not None:
            self.particle_scale = particle_scale
        if bond_scale is not None:
            self.bond_scale = bond_scale
        if scale is not None:
            self.scale = scale
        if color is not None:
            self.color = color
        if opacity is not None:
            self.opacity = opacity
        if edge_sharpness is not None:
            self.edge_sharpness = edge_sharpness

    @serialize_as("particle.scale")
    @properties.auto
    def particle_scale(self) -> float:
        """Scaling factor to apply to the balls."""

    @serialize_as("bond.scale")
    @properties.auto
    def bond_scale(self) -> float:
        """Scaling factor to apply to the bonds."""

    @serialize_as("scale")
    @properties.auto
    def scale(self) -> float:
        """"""

    @serialize_as("color")
    @properties.auto
    def color(self) -> SingleColor:
        """"""

    @serialize_as("opacity")
    @properties.auto
    def opacity(self) -> float:
        """"""

    @serialize_as("edge.sharpness")
    @properties.auto
    def edge_sharpness(self) -> float:
        """"""


class RenderCycles(Render):
    """Liquorice renderer with solid filled-in cycles."""

    _subgraph_ids: ClassVar[List[str]] = ["cycles", "cycle", "rings"]

    def __init__(
        self,
        *,
        scale: Optional[float] = None,
        cycle_scale: Optional[float] = None,
        color: Optional[SingleColor] = None,
        opacity: Optional[float] = None,
        bond_scale: Optional[float] = None,
        edge_sharpness: Optional[float] = None,
        **kwargs: Serializable
    ):
        """
        Liquorice renderer with solid filled-in cycles.

        :param scale:
        :param cycle_scale:
        :param color:
        :param opacity:
        :param bond_scale:
        :param edge_sharpness:
        """
        super().__init__(**kwargs)
        if scale is not None:
            self.scale = scale
        if cycle_scale is not None:
            self.cycle_scale = cycle_scale
        if color is not None:
            self.color = color
        if opacity is not None:
            self.opacity = opacity
        if bond_scale is not None:
            self.bond_scale = bond_scale
        if edge_sharpness is not None:
            self.edge_sharpness = edge_sharpness

    @serialize_as("scale")
    @properties.auto
    def scale(self) -> float:
        """"""

    @serialize_as("cycle.scale")
    @properties.auto
    def cycle_scale(self) -> float:
        """"""

    @serialize_as("color")
    @properties.auto
    def color(self) -> SingleColor:
        """"""

    @serialize_as("opacity")
    @properties.auto
    def opacity(self) -> float:
        """"""

    @serialize_as("bond.scale")
    @properties.auto
    def bond_scale(self) -> float:
        """"""

    @serialize_as("edge.sharpness")
    @properties.auto
    def edge_sharpness(self) -> float:
        """"""


class RenderRibbon(Render):
    """"""

    _subgraph_ids: ClassVar[List[str]] = ["elliptic spline", "ribbon", "ribbon spline"]

    def __init__(
        self,
        *,
        particle_widths: Optional[List[float]] = None,
        scale: Optional[float] = None,
        width: Optional[float] = None,
        speed: Optional[float] = None,
        input: Optional[List[SingleColor]] = None,
        color: Optional[SingleColor] = None,
        opacity: Optional[float] = None,
        edge_sharpness: Optional[float] = None,
        **kwargs: Serializable
    ):
        """


        :param particle_widths:
        :param scale:
        :param width:
        :param speed:
        :param input:
        :param color:
        :param opacity:
        :param edge_sharpness:
        """
        super().__init__(**kwargs)
        if particle_widths is not None:
            self.particle_widths = particle_widths
        if scale is not None:
            self.scale = scale
        if width is not None:
            self.width = width
        if speed is not None:
            self.speed = speed
        if input is not None:
            self.input = input
        if color is not None:
            self.color = color
        if opacity is not None:
            self.opacity = opacity
        if edge_sharpness is not None:
            self.edge_sharpness = edge_sharpness

    @serialize_as("particle.widths")
    @properties.auto
    def particle_widths(self) -> List[float]:
        """"""

    @serialize_as("scale")
    @properties.auto
    def scale(self) -> float:
        """"""

    @serialize_as("width")
    @properties.auto
    def width(self) -> float:
        """"""

    @serialize_as("speed")
    @properties.auto
    def speed(self) -> float:
        """"""

    @serialize_as("input")
    @properties.auto
    def input(self) -> List[SingleColor]:
        """"""

    @serialize_as("color")
    @properties.auto
    def color(self) -> SingleColor:
        """"""

    @serialize_as("opacity")
    @properties.auto
    def opacity(self) -> float:
        """"""

    @serialize_as("edge.sharpness")
    @properties.auto
    def edge_sharpness(self) -> float:
        """"""


class RenderGeometricRibbon(Render):
    """"""

    _subgraph_ids: ClassVar[List[str]] = [
        "tetrahedral spline",
        "geometric spline",
        "geometric",
    ]

    def __init__(
        self,
        *,
        particle_widths: Optional[List[float]] = None,
        scale: Optional[float] = None,
        width: Optional[float] = None,
        color: Optional[SingleColor] = None,
        opacity: Optional[float] = None,
        bond_scale: Optional[float] = None,
        edge_sharpness: Optional[float] = None,
        **kwargs: Serializable
    ):
        """


        :param particle_widths:
        :param scale:
        :param width:
        :param color:
        :param opacity:
        :param bond_scale:
        :param edge_sharpness:
        """
        super().__init__(**kwargs)
        if particle_widths is not None:
            self.particle_widths = particle_widths
        if scale is not None:
            self.scale = scale
        if width is not None:
            self.width = width
        if color is not None:
            self.color = color
        if opacity is not None:
            self.opacity = opacity
        if bond_scale is not None:
            self.bond_scale = bond_scale
        if edge_sharpness is not None:
            self.edge_sharpness = edge_sharpness

    @serialize_as("particle.widths")
    @properties.auto
    def particle_widths(self) -> List[float]:
        """"""

    @serialize_as("scale")
    @properties.auto
    def scale(self) -> float:
        """"""

    @serialize_as("width")
    @properties.auto
    def width(self) -> float:
        """"""

    @serialize_as("color")
    @properties.auto
    def color(self) -> SingleColor:
        """"""

    @serialize_as("opacity")
    @properties.auto
    def opacity(self) -> float:
        """"""

    @serialize_as("bond.scale")
    @properties.auto
    def bond_scale(self) -> float:
        """"""

    @serialize_as("edge.sharpness")
    @properties.auto
    def edge_sharpness(self) -> float:
        """"""


class RenderGoodsell(Render):
    """Goodsell style renderer which draws unshaded spheres with outlines between distant residues."""

    _subgraph_ids: ClassVar[List[str]] = ["goodsell", "goodsel"]

    def __init__(
        self,
        *,
        scale: Optional[float] = None,
        color: Optional[SingleColor] = None,
        opacity: Optional[float] = None,
        **kwargs: Serializable
    ):
        """
        Goodsell style renderer which draws unshaded spheres with outlines between distant residues.

        :param scale:
        :param color:
        :param opacity:
        """
        super().__init__(**kwargs)
        if scale is not None:
            self.scale = scale
        if color is not None:
            self.color = color
        if opacity is not None:
            self.opacity = opacity

    @serialize_as("scale")
    @properties.auto
    def scale(self) -> float:
        """"""

    @serialize_as("color")
    @properties.auto
    def color(self) -> SingleColor:
        """"""

    @serialize_as("opacity")
    @properties.auto
    def opacity(self) -> float:
        """"""


class RenderHydrogenBondCaps(Render):
    """Render cones to indicate directions which hydrogen bonds should exist."""

    _subgraph_ids: ClassVar[List[str]] = ["hydrogen cap"]

    def __init__(
        self,
        *,
        acceptor_scale: Optional[float] = None,
        donor_scale: Optional[float] = None,
        acceptor_color: Optional[SingleColor] = None,
        donor_color: Optional[SingleColor] = None,
        scale: Optional[float] = None,
        color: Optional[SingleColor] = None,
        opacity: Optional[float] = None,
        **kwargs: Serializable
    ):
        """
        Render cones to indicate directions which hydrogen bonds should exist.

        :param acceptor_scale:
        :param donor_scale:
        :param acceptor_color:
        :param donor_color:
        :param scale:
        :param color:
        :param opacity:
        """
        super().__init__(**kwargs)
        if acceptor_scale is not None:
            self.acceptor_scale = acceptor_scale
        if donor_scale is not None:
            self.donor_scale = donor_scale
        if acceptor_color is not None:
            self.acceptor_color = acceptor_color
        if donor_color is not None:
            self.donor_color = donor_color
        if scale is not None:
            self.scale = scale
        if color is not None:
            self.color = color
        if opacity is not None:
            self.opacity = opacity

    @serialize_as("acceptor.scale")
    @properties.auto
    def acceptor_scale(self) -> float:
        """"""

    @serialize_as("donor.scale")
    @properties.auto
    def donor_scale(self) -> float:
        """"""

    @serialize_as("acceptor.color")
    @properties.auto
    def acceptor_color(self) -> SingleColor:
        """"""

    @serialize_as("donor.color")
    @properties.auto
    def donor_color(self) -> SingleColor:
        """"""

    @serialize_as("scale")
    @properties.auto
    def scale(self) -> float:
        """"""

    @serialize_as("color")
    @properties.auto
    def color(self) -> SingleColor:
        """"""

    @serialize_as("opacity")
    @properties.auto
    def opacity(self) -> float:
        """"""


class RenderHyperballs(Render):
    """Draws spheres connected by hyperboloids."""

    _subgraph_ids: ClassVar[List[str]] = ["hyperballs", "hyperball"]

    def __init__(
        self,
        *,
        scale: Optional[float] = None,
        tension: Optional[float] = None,
        color: Optional[SingleColor] = None,
        opacity: Optional[float] = None,
        bond_scale: Optional[float] = None,
        edge_sharpness: Optional[float] = None,
        **kwargs: Serializable
    ):
        """
        Draws spheres connected by hyperboloids.

        :param scale:
        :param tension:
        :param color:
        :param opacity:
        :param bond_scale:
        :param edge_sharpness:
        """
        super().__init__(**kwargs)
        if scale is not None:
            self.scale = scale
        if tension is not None:
            self.tension = tension
        if color is not None:
            self.color = color
        if opacity is not None:
            self.opacity = opacity
        if bond_scale is not None:
            self.bond_scale = bond_scale
        if edge_sharpness is not None:
            self.edge_sharpness = edge_sharpness

    @serialize_as("scale")
    @properties.auto
    def scale(self) -> float:
        """"""

    @serialize_as("tension")
    @properties.auto
    def tension(self) -> float:
        """"""

    @serialize_as("color")
    @properties.auto
    def color(self) -> SingleColor:
        """"""

    @serialize_as("opacity")
    @properties.auto
    def opacity(self) -> float:
        """"""

    @serialize_as("bond.scale")
    @properties.auto
    def bond_scale(self) -> float:
        """"""

    @serialize_as("edge.sharpness")
    @properties.auto
    def edge_sharpness(self) -> float:
        """"""


class RenderLiquorice(Render):
    """Render stick representation with rounded bonds."""

    _subgraph_ids: ClassVar[List[str]] = ["liquorice", "licorice"]

    def __init__(
        self,
        *,
        scale: Optional[float] = None,
        color: Optional[SingleColor] = None,
        opacity: Optional[float] = None,
        bond_scale: Optional[float] = None,
        edge_sharpness: Optional[float] = None,
        **kwargs: Serializable
    ):
        """
        Render stick representation with rounded bonds.

        :param scale:
        :param color:
        :param opacity:
        :param bond_scale:
        :param edge_sharpness:
        """
        super().__init__(**kwargs)
        if scale is not None:
            self.scale = scale
        if color is not None:
            self.color = color
        if opacity is not None:
            self.opacity = opacity
        if bond_scale is not None:
            self.bond_scale = bond_scale
        if edge_sharpness is not None:
            self.edge_sharpness = edge_sharpness

    @serialize_as("scale")
    @properties.auto
    def scale(self) -> float:
        """"""

    @serialize_as("color")
    @properties.auto
    def color(self) -> SingleColor:
        """"""

    @serialize_as("opacity")
    @properties.auto
    def opacity(self) -> float:
        """"""

    @serialize_as("bond.scale")
    @properties.auto
    def bond_scale(self) -> float:
        """"""

    @serialize_as("edge.sharpness")
    @properties.auto
    def edge_sharpness(self) -> float:
        """"""


class RenderNoodles(Render):
    """Liquorice-like renderer with bonds curved based on neighbouring bonds."""

    _subgraph_ids: ClassVar[List[str]] = ["noodles"]

    def __init__(
        self,
        *,
        scale: Optional[float] = None,
        color: Optional[SingleColor] = None,
        opacity: Optional[float] = None,
        edge_sharpness: Optional[float] = None,
        **kwargs: Serializable
    ):
        """
        Liquorice-like renderer with bonds curved based on neighbouring bonds.

        :param scale:
        :param color:
        :param opacity:
        :param edge_sharpness:
        """
        super().__init__(**kwargs)
        if scale is not None:
            self.scale = scale
        if color is not None:
            self.color = color
        if opacity is not None:
            self.opacity = opacity
        if edge_sharpness is not None:
            self.edge_sharpness = edge_sharpness

    @serialize_as("scale")
    @properties.auto
    def scale(self) -> float:
        """"""

    @serialize_as("color")
    @properties.auto
    def color(self) -> SingleColor:
        """"""

    @serialize_as("opacity")
    @properties.auto
    def opacity(self) -> float:
        """"""

    @serialize_as("edge.sharpness")
    @properties.auto
    def edge_sharpness(self) -> float:
        """"""


class RenderPeptidePlanes(Render):
    """Render solid planes aligned with the peptide bonds."""

    _subgraph_ids: ClassVar[List[str]] = ["peptide planes"]

    def __init__(
        self,
        *,
        scale: Optional[float] = None,
        color: Optional[SingleColor] = None,
        opacity: Optional[float] = None,
        bond_scale: Optional[float] = None,
        edge_sharpness: Optional[float] = None,
        **kwargs: Serializable
    ):
        """
        Render solid planes aligned with the peptide bonds.

        :param scale:
        :param color:
        :param opacity:
        :param bond_scale:
        :param edge_sharpness:
        """
        super().__init__(**kwargs)
        if scale is not None:
            self.scale = scale
        if color is not None:
            self.color = color
        if opacity is not None:
            self.opacity = opacity
        if bond_scale is not None:
            self.bond_scale = bond_scale
        if edge_sharpness is not None:
            self.edge_sharpness = edge_sharpness

    @serialize_as("scale")
    @properties.auto
    def scale(self) -> float:
        """"""

    @serialize_as("color")
    @properties.auto
    def color(self) -> SingleColor:
        """"""

    @serialize_as("opacity")
    @properties.auto
    def opacity(self) -> float:
        """"""

    @serialize_as("bond.scale")
    @properties.auto
    def bond_scale(self) -> float:
        """"""

    @serialize_as("edge.sharpness")
    @properties.auto
    def edge_sharpness(self) -> float:
        """"""


class RenderSpheres(Render):
    """Render spheres without bonds."""

    _subgraph_ids: ClassVar[List[str]] = [
        "spheres",
        "sphere",
        "ball",
        "balls",
        "space filling",
    ]

    def __init__(
        self,
        *,
        particle_scale: Optional[float] = None,
        outline_depth: Optional[float] = None,
        outline_width: Optional[float] = None,
        scale: Optional[float] = None,
        color: Optional[SingleColor] = None,
        opacity: Optional[float] = None,
        **kwargs: Serializable
    ):
        """
        Render spheres without bonds.

        :param particle_scale:
        :param outline_depth:
        :param outline_width:
        :param scale:
        :param color:
        :param opacity:
        """
        super().__init__(**kwargs)
        if particle_scale is not None:
            self.particle_scale = particle_scale
        if outline_depth is not None:
            self.outline_depth = outline_depth
        if outline_width is not None:
            self.outline_width = outline_width
        if scale is not None:
            self.scale = scale
        if color is not None:
            self.color = color
        if opacity is not None:
            self.opacity = opacity

    @serialize_as("particle.scale")
    @properties.auto
    def particle_scale(self) -> float:
        """"""

    @serialize_as("outline.depth")
    @properties.auto
    def outline_depth(self) -> float:
        """"""

    @serialize_as("outline.width")
    @properties.auto
    def outline_width(self) -> float:
        """"""

    @serialize_as("scale")
    @properties.auto
    def scale(self) -> float:
        """"""

    @serialize_as("color")
    @properties.auto
    def color(self) -> SingleColor:
        """"""

    @serialize_as("opacity")
    @properties.auto
    def opacity(self) -> float:
        """"""


class RenderTube(Render):
    """Render solid curved tube through alpha carbons."""

    _subgraph_ids: ClassVar[List[str]] = [
        "spline",
        "circular spline",
        "cylindrical spline",
    ]

    def __init__(
        self,
        *,
        scale: Optional[float] = None,
        color: Optional[SingleColor] = None,
        opacity: Optional[float] = None,
        edge_sharpness: Optional[float] = None,
        **kwargs: Serializable
    ):
        """
        Render solid curved tube through alpha carbons.

        :param scale:
        :param color:
        :param opacity:
        :param edge_sharpness:
        """
        super().__init__(**kwargs)
        if scale is not None:
            self.scale = scale
        if color is not None:
            self.color = color
        if opacity is not None:
            self.opacity = opacity
        if edge_sharpness is not None:
            self.edge_sharpness = edge_sharpness

    @serialize_as("scale")
    @properties.auto
    def scale(self) -> float:
        """"""

    @serialize_as("color")
    @properties.auto
    def color(self) -> SingleColor:
        """"""

    @serialize_as("opacity")
    @properties.auto
    def opacity(self) -> float:
        """"""

    @serialize_as("edge.sharpness")
    @properties.auto
    def edge_sharpness(self) -> float:
        """"""
