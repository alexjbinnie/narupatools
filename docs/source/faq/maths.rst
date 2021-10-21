***************
Maths & Physics
***************

Narupatools contains a suite of methods for performing physical calculations.

How do I...
-----------

... easily create a vector?
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :obj:`narupatools.physics.vector.vector` method creates a one-dimensional numpy array with the float type:

.. testcode:: vector

   from narupatools.physics.vector import vector

   vec = vector(0, 1, 3)  # NumPy array of dtype float

   vec = np.array([0, 1, 3])  # NumPy array of dtype int (may not be desired)

   vec = np.array([0, 1, 3], dtype=float)  # NumPy array of dtype float (vector() is shorthand for this)

... take dot products of two vectors (or sets of vectors)?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :obj:`narupatools.physics.vector.dot_product` method takes the dot product of two vectors. It handles higher dimensions than the builtin :obj:`numpy.dot` function, which only works for up to 2D vectors.

... take cross products of two vectors (or arrays of vectors)?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :obj:`narupatools.physics.vector.cross_product` method takes the cross product of two vectors.

... get the magnitude of a vector (or array of vectors)?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :obj:`narupatools.physics.vector.magnitude` function gets the norm of a vector or set of vectors.

.. testcode:: dynamics_run_steps

   # Run forever in the background. Set block=True to run in the current thread.
   dynamics.run(block=False)

.. testcleanup:: dynamics_run_steps

   dynamics.stop()

... normalize vectors or quaternions?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :obj:`narupatools.physics.vector.normalized` function can normalize vectors, quaternions and arrays of vectors or quaternions.

... project one vector onto another?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The vector projection of a vector :math:`a` onto a vector :math:`b` is the component of :math:`a` which is parallel to :math:`b`. This can be obtained using the :obj:`narupatools.physics.vector.vector_projection` method.

... reject one vector from another?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The vector rejection of a vector :math:`a` onto a vector :math:`b` is the component of :math:`a` which is perpendicular to :math:`b`. This can be obtained using the :obj:`narupatools.physics.vector.vector_rejection` method.

... get the distance between two points?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The distance between two points can be obtained using :obj:`narupatools.physics.vector.distance`. If the square of the distance is needed, use :obj:`narupatools.physics.vector.sqr_distance` to avoid calculating square roots.

... get the angle between two vectors?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The angle in radians between two vectors can be obtained using :obj:`narupatools.physics.vector.angle`.