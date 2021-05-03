################
HDF5 File Format
################

HDF5 is a file format for storing large amounts of data. Like XML, HDF5 merely describes a format for storing information such as arrays. Hence, to define a format to store trajectories in HDF5, we must define exactly how our data will be stored.

The following format is an extension of the MDTraj trajectory format, as defined `here <https://mdtraj.org/1.9.4/hdf5_format.html>`_. See that article for the reasoning behind the choice of format.

NarupaTools Trajectory 1.0
==========================

The NarupaTools trajectory format is a superset of the MDTraj trajectory 1.1 convention.

Global Attributes
-----------------

* **Conventions (required)**: Space/comma separated list of tokens defining conventions followed by the file. MDTraj requires that 'Pande' be one of the conventions present. NarupaTools additionally requires that 'NarupaTools' also be present to indicate a trajectory adheres to this standard.

* **ConventionVersion (required)**: Should be set to "1.1" to indicate the convention version for MDTraj.

* **NarupaToolsConventionVersion (required)**: Should be set to "1.0", to indicate the convention version for NarupaTools.

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
