from openmm import CustomIntegrator

class MTSIntegrator(CustomIntegrator):
    def __init__(self, dt: float, groups): ...

class MTSLangevinIntegrator(CustomIntegrator):
    def __init__(self, temperature: float, friction: float, dt: float, groups): ...
