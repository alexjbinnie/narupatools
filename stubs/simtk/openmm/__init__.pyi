import logging

from openmm import *

from . import app

logging.getLogger().warning(
    "Warning: importing 'simtk.openmm' is deprecated.  Import 'openmm' instead."
)
