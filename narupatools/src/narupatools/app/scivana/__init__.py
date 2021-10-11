# This file is part of narupatools (https://github.com/alexjbinnie/narupatools).
# Copyright (c) Alex Jamieson-Binnie. All rights reserved.
#
# narupatools is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# narupatools is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with narupatools.  If not, see <http://www.gnu.org/licenses/>.

"""Code relating to using Scivana."""

from ._autogen._color import (
    Color,
    ColorByElement,
    ColorByFloatGradient,
    ColorByParticleIndex,
    ColorByParticleType,
    ColorByResidueIndex,
    ColorByResidueName,
    ColorBySecondaryStructure,
)
from ._autogen._render import (
    Render,
    RenderBallAndStick,
    RenderCycles,
    RenderGeometricRibbon,
    RenderGoodsell,
    RenderHydrogenBondCaps,
    RenderHyperballs,
    RenderLiquorice,
    RenderNoodles,
    RenderPeptidePlanes,
    RenderRibbon,
    RenderSpheres,
    RenderTube,
)
from ._autogen._renderer import Renderer
from ._autogen._scale import (
    Scale,
    ScaleByCovalent,
    ScaleBySecondaryStructure,
    ScaleByVanDerWaals,
)
from ._autogen._sequence import Sequence, SequenceByPolypeptideAlphaCarbons
from ._camera import CameraView
from ._particle_selection import ParticleSelection
from ._particle_visualiser import ParticleVisualisation

__all__ = [
    "Color",
    "ColorByElement",
    "ColorByFloatGradient",
    "ColorByParticleIndex",
    "ColorByParticleType",
    "ColorByResidueName",
    "ColorByResidueIndex",
    "ColorBySecondaryStructure",
    "Render",
    "RenderBallAndStick",
    "RenderHyperballs",
    "Renderer",
    "RenderTube",
    "RenderRibbon",
    "RenderLiquorice",
    "RenderCycles",
    "RenderGoodsell",
    "RenderNoodles",
    "RenderSpheres",
    "RenderPeptidePlanes",
    "RenderHydrogenBondCaps",
    "RenderGeometricRibbon",
    "Scale",
    "Sequence",
    "SequenceByPolypeptideAlphaCarbons",
    "ScaleByCovalent",
    "ScaleByVanDerWaals",
    "ScaleBySecondaryStructure",
    "ParticleVisualisation",
    "ParticleSelection",
    "CameraView",
]
