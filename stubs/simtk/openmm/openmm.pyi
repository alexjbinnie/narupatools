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

from typing import Sequence, Union, overload, List

import numpy as np
from simtk.openmm.vec3 import Vec3
from simtk.unit import Quantity
from typing_extensions import Literal

Vector3Array = Union[
    np.ndarray, Quantity[np.ndarray], Sequence[Vec3], Quantity[Sequence[Vec3]]
]

class Context:
    def getState(
        self,
        getPositions: bool = ...,
        getForces: bool = ...,
        getEnergy: bool = ...,
        getVelocities: bool = ...,
        getParameters: bool = ...,
    ) -> State: ...
    def setPositions(self, positions: Vector3Array) -> None: ...
    def setVelocities(self, positions: Vector3Array) -> None: ...
    def getSystem(self) -> System: ...
    def reinitialize(self, preserveState: bool = ...) -> None: ...

class Integrator:
    def getStepSize(self) -> Quantity[float]: ...
    def setStepSize(self, float: Union[float, Quantity[float]]) -> None: ...
    def step(self, steps: int) -> None: ...

class CustomIntegrator(Integrator):
    def __init__(self, stepSize: float): ...
    def addPerDofVariable(self, name: str, initialValue: float) -> None: ...
    def addUpdateContextState(self) -> None: ...
    def addComputePerDof(self, variable: str, expression: str) -> None: ...
    def addConstrainPositions(self) -> None: ...
    def addConstrainVelocities(self) -> None: ...

class State:
    @overload
    def getForces(self, asNumpy: Literal[True]) -> Quantity[np.ndarray]: ...
    @overload
    def getForces(self, asNumpy: Literal[False]) -> Quantity[Sequence[Vec3]]: ...
    @overload
    def getForces(
        self, asNumpy: bool = ...
    ) -> Union[Quantity[np.ndarray], Quantity[Sequence[Vec3]]]: ...
    def getPotentialEnergy(self) -> Quantity[float]: ...
    def getKineticEnergy(self) -> Quantity[float]: ...
    @overload
    def getPositions(self, asNumpy: Literal[True]) -> Quantity[np.ndarray]: ...
    @overload
    def getPositions(self, asNumpy: Literal[False]) -> Quantity[Sequence[Vec3]]: ...
    @overload
    def getPositions(
        self, asNumpy: bool = ...
    ) -> Union[Quantity[np.ndarray], Quantity[Sequence[Vec3]]]: ...
    @overload
    def getVelocities(self, asNumpy: Literal[True]) -> Quantity[np.ndarray]: ...
    @overload
    def getVelocities(self, asNumpy: Literal[False]) -> Quantity[Sequence[Vec3]]: ...
    @overload
    def getVelocities(
        self, asNumpy: bool = ...
    ) -> Union[Quantity[np.ndarray], Quantity[Sequence[Vec3]]]: ...
    def getPeriodicBoxVectors(self) -> Quantity[Sequence[Vec3]]: ...

class System:
    def getNumParticles(self) -> int: ...
    def addParticle(self, mass: float) -> int: ...
    def getParticleMass(self, index: int) -> Quantity[float]: ...
    def setParticleMass(self, index: int, mass: float) -> None: ...
    def getNumConstraints(self) -> int: ...
    def usesPeriodicBoundaryConditions(self) -> bool: ...
    def getDefaultPeriodicBoxVectors(self) -> Quantity[Sequence[Vec3]]: ...
    def addForce(self, force: Force) -> int: ...
    def getForces(self) -> List[Force]: ...

class LangevinIntegrator(Integrator):
    def __init__(
        self,
        temperature: Union[float, Quantity[float]],
        frictionCoeff: Union[float, Quantity[float]],
        stepSize: Union[float, Quantity[float]],
    ): ...
    def getFriction(self) -> Quantity[float]: ...

class Force: ...
class TabulatedFunction: ...

class CustomExternalForce(Force):
    def __init__(self, energy : str): ...
    def setParticleParameters(
        self, index: int, particle: int, parameters: Union[Sequence[float], np.ndarray]
    ) -> None: ...
    def updateParametersInContext(self, context: Context) -> None: ...
    def addPerParticleParameter(self, name: str) -> int: ...
    def addParticle(self, particle: int, parameters: Sequence[float]): ...
    def getEnergyFunction(self) -> str: ...
    def getNumParticles(self) -> int: ...

class XmlSerializer:
    @staticmethod
    def serialize(
        object: Union[System, Force, Integrator, State, TabulatedFunction]
    ) -> str: ...
    @staticmethod
    def deserialize(
        inputString: str,
    ) -> Union[System, Force, Integrator, State, TabulatedFunction]: ...
