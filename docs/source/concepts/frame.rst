###################
Narupa Frame Format
###################

Particle Positions
------------------

Key
    `particle.positions`

Canonical representation
    Array of dimensions :math:`(N_{\text{particles}}, 3)`, storing the x, y and z coordinates of each particle.

Dimension
    `length`

Standard Units
    `nanometers`

Python Representation
    NumPy float array with shape `(N, 3)`.

FrameData Representation
    Float array of dimension `3N`, storing the coordinates in the form :math:`(x_0, y_0, z_0, x_1, y_1, z_1, ...)`, under the key `particle.positions`.


Particle Velocities
-------------------

Key
    `particle.velocities`

Canonical representation
    Array of dimensions :math:`(N_{\text{particles}}, 3)`, storing the x, y and z components of the velocity of each particle.

Dimension
    `length / time`

Standard Units
    `nanometers / picoseconds`

Python Representation
    NumPy float array with shape `(N, 3)`.

FrameData Representation
    Float array of dimension `3N`, storing the components in the form :math:`(v_0^x, v_0^y, v_0^z, v_1^x, v_1^y, v_1^z, ...)`, under the key `particle.velocities`.


Particle Forces
---------------

Key
    `particle.forces`

Canonical representation
    Array of dimensions :math:`(N_{\text{particles}}, 3)`, storing the x, y and z components of the force on each particle.

Dimension
    `length * mass / (time * time)`

Standard Units
    `kilojoules / mole / picoseconds`

Python Representation
    NumPy float array with shape `(N, 3)`.

FrameData Representation
    Float array of dimension `3N`, storing the components in the form :math:`(F_0^x, F_0^y, F_0^z, F_1^x, F_1^y, F_1^z, ...)`, under the key `particle.forces`.


Particle Element (Atomic Number)
--------------------------------

Key
    `particle.elements`

Canonical representation
    Array of dimensions :math:`(N_{\text{particles}})`, storing the atomic number `Z` of each particle. A value of `0` shall be taken to indicate a particle which does not have an element.

Dimension
    N/A

Standard Units
    N/A

Python Representation
    NumPy int array with shape `(N,)`.

FrameData Representation
    Index array of dimension `N`, under the key `particle.elements`.


Particle Residues
-----------------

Key
    `particle.residues`

Canonical representation
    Array of dimensions :math:`(N_{\text{particles}})`, storing the index of the residue this particle belongs to.

Dimension
    N/A

Standard Units
    N/A

Python Representation
    NumPy int array with shape `(N,)`.

FrameData Representation
    Index array of dimension `N`, under the key `particle.residues`.


Particle Names
--------------

Key
    `particle.names`

Canonical representation
    Array of dimensions :math:`(N_{\text{particles}})`, storing a string name for each particle. The name has no meaning to the protocol. It is usually defined by the force field.

Dimension
    N/A

Standard Units
    N/A

Python Representation
    NumPy object array with shape `(N,)`.

FrameData Representation
    String array of dimension `N`, under the key `particle.names`.


Particle Types
--------------

Key
    `particle.types`

Canonical representation
    Array of dimensions :math:`(N_{\text{particles}})`, storing a string type for each particle. The type has no meaning to the protocol. Its use is for identifying particles when they do not correspond to a specific element, for example when running coarse grained simulations.

Dimension
    N/A

Standard Units
    N/A

Python Representation
    NumPy object array with shape `(N,)`.

FrameData Representation
    String array of dimension `N`, under the key `particle.types`.


Particle Masses
---------------

Key
    `particle.masses`

Canonical representation
    Array of dimensions :math:`(N_{\text{particles}})`, storing the mass `m` of each particle.

Dimension
    `mass`

Standard Units
    `dalton`

Python Representation
    NumPy float array with shape `(N,)`.

FrameData Representation
    Float array of dimension `N`, under the key `particle.masses`.

Calculation
    Implementations may allow the value of this key to be inferred from the existence of the `particle.elements` key. If this is present, masses may be derived from the masses of each particle.


Particle Charges
----------------

Key
    `particle.charges`

Canonical representation
    Array of dimensions :math:`(N_{\text{particles}})`, storing the charge `q` of each particle.

Dimension
    `charge`

Standard Units
    `elementary charge`

Python Representation
    NumPy float array with shape `(N,)`.

FrameData Representation
    Float array of dimension `N`, under the key `particle.charges`.


Particle Count
----------------

Key
    `particle.count`

Canonical representation
    Integer :math:`N_{\text{particles}}`, denoting the number of particles.

Dimension
    N/A

Standard Units
    N/A

Python Representation
    `int`

FrameData Representation
    NumberValue, under the key `particle.count`.


Residue Names
-------------

Key
    `residue.names`

Canonical representation
    Array of dimensions :math:`(N_{\text{residues}})`, storing a string name for each residue.

Dimension
    N/A

Standard Units
    N/A

Python Representation
    NumPy object array with shape `(N,)`.

FrameData Representation
    String array of dimension `N`, under the key `residue.names`.


Residue IDs
-----------

Key
    `residue.ids`

Canonical representation
    Array of dimensions :math:`(N_{\text{residues}})`, storing an id for each residue. This is often known as the seg id.

Dimension
    N/A

Standard Units
    N/A

Python Representation
    NumPy object array with shape `(N,)`.

FrameData Representation
    String array of dimension `N`, under the key `residue.ids`.


Residue Chains
--------------

Key
    `residue.chains`

Canonical representation
    Array of length :math:`N_{\text{residues}}`, storing the index of the chain each residue belongs to.

Dimension
    N/A

Standard Units
    N/A

Python Representation
    NumPy int array with shape `(N,)`.

FrameData Representation
    Index array of dimension `N`, under the key `residue.chains`.


Residue Count
-------------

Key
    `residue.count`

Canonical representation
    Integer :math:`N_{\text{residues}}`, denoting the number of residues.

Dimension
    N/A

Standard Units
    N/A

Python Representation
    `int`

FrameData Representation
    NumberValue, under the key `residue.count`.


Chain Names
-----------

Key
    `chain.names`

Canonical representation
    Array of dimensions :math:`(N_{\text{chains}})`, storing a string name for each chain.

Dimension
    N/A

Standard Units
    N/A

Python Representation
    NumPy object array with shape `(N,)`.

FrameData Representation
    String array of dimension `N`, under the key `chain.names`.


Chain Count
-----------

Key
    `chain.count`

Canonical representation
    Integer :math:`N_{\text{chains}}`, denoting the number of residues.

Dimension
    N/A

Standard Units
    N/A

Python Representation
    `int`

FrameData Representation
    NumberValue, under the key `chain.count`.


Bond Pairs
----------

Key
    `bond.pairs`

Canonical representation
    Array of dimensions :math:`(N_{\text{bonds}}, 2)`, storing the indices of the two particles involved in each bond.

Dimension
    N/A

Standard Units
    N/A

Python Representation
    NumPy int array with shape `(N, 2)`.

FrameData Representation
    Index array of dimension `2N`, storing the indices in the form :math:`(a_0, b_0, a_1, b_1, ...)`, under the key `bond.pairs`.


Bond Orders
-----------

Key
    `bond.orders`

Canonical representation
    Array of length :math:`N_{\text{bonds}}`, storing the bond order of each bond.

Dimension
    N/A

Standard Units
    N/A

Python Representation
    NumPy int array with shape `(N,)`.

FrameData Representation
    Index array of dimension `N`, under the key `bond.orders`.


Bond Count
----------

Key
    `bond.count`

Canonical representation
    Integer :math:`N_{\text{bonds}}`, denoting the number of bonds.

Dimension
    N/A

Standard Units
    N/A

Python Representation
    `int`

FrameData Representation
    NumberValue, under the key `bond.count`.


Box Vectors
-----------

Key
    `box.vectors`

Canonical representation
    :math:`3\times 3` matrix, whose columns denote the three axes of the unit scell.

Dimension
    `length`

Standard Units
    `nanometer`

Python Representation
    NumPy float array of shape `(3,3)`.

FrameData Representation
    Float array of length `9`, with components stored in the order :math:`(a_x, a_y, a_z, b_x, b_y, b_z, c_x, c_y, c_z)`.


Potential Energy
----------------

Key
    `energy.potential`

Canonical representation
    Decimal :math:`PE`, denoting the total potential energy of the system.

Dimension
    `mass * length * length / (time * time)`

Standard Units
    `kilojoule / mole`

Python Representation
    `float`

FrameData Representation
    NumberValue, under the key `energy.potential`.


Kinetic Energy
--------------

Key
    `energy.kinetic`

Canonical representation
    Decimal :math:`KE`, denoting the total kinetic energy of the system.

Dimension
    `mass * length * length / (time * time)`

Standard Units
    `kilojoule / mole`

Python Representation
    `float`

FrameData Representation
    NumberValue, under the key `energy.kinetic`.

Calculation
    Kinetic energy may be calculated if both `particle.masses` and `particle.velocities` are present. The kinetic energy is given by:

    .. math:: KE = \sum_i \frac{1}{2} m_i v_i^2


Simulation Elapsed Time
-----------------------

Key
    `simulation.elapsed_time`

Canonical representation
    Decimal :math:`t`, denoting the time that has passed in simulation time. This time resets when the simulation is reset.

Dimension
    `time`

Standard Units
    `picosecond`

Python Representation
    `float`

FrameData Representation
    NumberValue, under the key `simulation.elapsed_time`.


Simulation Total Time
---------------------

Key
    `simulation.total_time`

Canonical representation
    Decimal :math:`t`, denoting the time that has passed in simulation time. This time does not reset when the simulation is reset.

Dimension
    `time`

Standard Units
    `picosecond`

Python Representation
    `float`

FrameData Representation
    NumberValue, under the key `simulation.total_time`.


Simulation Elapsed Steps
------------------------

Key
    `simulation.elapsed_steps`

Canonical representation
    Integer :math:`n`, denoting the number of steps that has passed in the simulation. This resets when the simulation is reset.

Dimension
    N/A

Standard Units
    N/A

Python Representation
    `int`

FrameData Representation
    NumberValue, under the key `simulation.elapsed_steps`.


Simulation Total Steps
----------------------

Key
    `simulation.total_steps`

Canonical representation
    Integer :math:`n`, denoting the number of steps that has passed in the simulation. This does not reset when the simulation is reset.

Dimension
    N/A

Standard Units
    N/A

Python Representation
    `int`

FrameData Representation
    NumberValue, under the key `simulation.total_steps`.


Derived Particle Properties
===========================

The following are fields which should not be stored to disk, but they could be computed from other fields.


Particle Momenta
----------------

Key
    `particle.momenta`

Canonical representation
    Array of dimensions :math:`(N_{\text{particles}})`, storing the x, y and z components of the momentum :math:`p` of each particle.

Dimension
    `mass * length / time`

Standard Units
    `dalton * nanometer / picosecond`

Calculation
    Momenta may be calculated if both `particle.masses` and `particle.velocities` are present. The momenta may be calculated for each particle by:

    .. math:: p_i = m_i v_i


Particle Acceleration
---------------------

Key
    `particle.accelerations`

Canonical representation
    Array of dimensions :math:`(N_{\text{particles}}, 3)`, storing the x, y and z components of the acceleration :math:`a` of each particle.

Dimension
    `length / (time * time)`

Standard Units
    `nanometer / (picosecond * picosecond)`

Calculation
    Acceleration may be calculated if both `particle.masses` and `particle.forces` are present. The acceleration may be calculated for each particle by:

    .. math:: a_i = F_i / m_i