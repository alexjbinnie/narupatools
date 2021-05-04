################
HDF5 File Format
################

HDF5 is a file format for storing large amounts of data. Like XML, HDF5 merely describes a format for storing information such as arrays. Hence, to define a format to store trajectories in HDF5, we must define exactly how our data will be stored.

The following format is an extension of the MDTraj trajectory format, as defined `here <https://mdtraj.org/1.9.4/hdf5_format.html>`_. See that article for the reasoning behind the choice of format.

NarupaTools Trajectory 1.0
==========================

The NarupaTools 1.0 trajectory format is a superset of the MDTraj trajectory 1.1 convention.

Global Attributes
-----------------

* **conventions (required)**: Space/comma separated list of tokens defining conventions followed by the file. MDTraj requires that 'Pande' be one of the conventions present. NarupaTools additionally requires that 'NarupaTools' also be present to indicate a trajectory adheres to this standard.

* **conventionVersion (required)**: Should be set to "1.1" to indicate the convention version for MDTraj.

* **narupaToolsConventionVersion (required)**: Should be set to "1.0", to indicate the convention version for *narupatools*.

* **program (required)**: Name of the creating program or module.

* **programVersion (required)**: Version of the creating program or module.

Arrays
------

Within the root group '/' the following arrays can be specified. Each array should have an attribute with the key 'units' with a string value specifying the units.

* **coordinates (required)**: Cartesian coordinates of the particles.

    * Shape: (n_frames, n_atoms, 3)
    * Type: float32
    * Standard Units: ”nanometers”

* **time (optional)**: Time of each frame.

    * Shape: (n_frames)
    * Type: float32
    * Standard Units: "picoseconds"

* **cell_lengths (optional)**: Cell lengths of the unit cell. The cell is standardised, with :math:`a` lying on the x axis and :math:`b` in the x-y plane. The origin of the unit cell lies at :math:`(0,0,0)`. If any dimension is not periodic, then the cell length should be set to 0.

    * Shape: (n_frames, 3)
    * Type: float32
    * Standard Units: "nanometers"

* **cell_angles (optional)**: Cell angles of the unit cell.

    * Shape: (n_frames, 3)
    * Type: float32
    * Standard Units: "degrees"

* **velocities (optional)**: Velocities of the particles.

    * Shape: (n_frames, n_atoms, 3)
    * Type: float32
    * Standard Units: ”nanometers/picosecond”

* **forces (optional)**: Forces on each particle.

    * Shape: (n_frames, n_atoms, 3)
    * Type: float32
    * Standard Units: "kJ/mol/nanometer"

* **kineticEnergy (optional)**: Kinetic energy of each frame.

    * Shape: (n_frames)
    * Type: float32
    * Standard Units: "kJ/mol"

* **potentialEnergy (optional)**: Potential energy of each frame.

    * Shape: (n_frames)
    * Type: float32
    * Standard Units: "kJ/mol"

* **topology (optional)**: Topology specified in ASCII encoded JSON. See the MDTraj specification for its format.

    * Shape: (1, ...)
    * Type: string

Interactions
============

The *narupatools* HDF5 format additionally stores interactions applied to the system. They are stored under a group 'interactions' with the key '/interactions'.

Each interaction is logged as a separate group under this group. For example, an interaction could be saved under '/interactions/interaction-my-id'. The key that each interaction is stored under is identical to the key that it was stored in the shared state.

Interaction Attributes
----------------------

Each interaction has the following attributes:

* **type (required)**: Type of the interaction, such as 'spring' or 'gaussian'.

* **startIndex (required)**: The index of the first frame this interaction was applied to.

* **endIndex (required)**: The index of the last frame this interaction was applied to.

Interaction Arrays
------------------

* **indices (required)**: Indices of the particles involved in the interaction.

    * Shape: (n_atoms_interaction)
    * Type: int32
    * Standard Units: "kJ/mol"

* **position (required)**: Position of the interaction.

    * Shape: (n_frames, 3)
    * Type: float32
    * Standard Units: "nanometers"

* **forces (required)**: Force on each particle due to the interaction at each frame.

    * Shape: (n_frames, n_atoms_interactions, 3)
    * Type: float32
    * Standard Units: "kJ/mol/nanometer"

* **potentialEnergy (required)**: Potential energy of the interaction at each frame.

    * Shape: (n_frames)
    * Type: float32
    * Standard Units: "kJ/mol"

* **frameIndex (required)**: Index of the corresponding frame in the main trajectory.

    * Shape: (n_frames)
    * Type: int32

* **scale (required)**: Interaction scale at each frame

    * Shape: (n_frames)
    * Type: float32