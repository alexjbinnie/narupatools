##############################
Interactive Molecular Dynamics
##############################

Narupa's main feature is that it enabled interactive molecular dynamics to be performed on simulations running in various packages.

By an **interaction**, we mean a force applied to one or more atoms in a simulations, normally for a short period of time. Interactions are applied to a group of atoms by 'mass-weighting'. The group of atoms are treated as a single composite point particle, with total mass :math:`M` and center of mass :math:`r_{\text{COM}}`:

.. math:: M = \sum_i m_i

.. math:: r_{\text{COM}} = \frac{1}{M} \sum m_i r_i

The force :math:`F` to apply to this 'composite particle' is calculated depending on the type of interaction. Then this force is distributed across each of the affected particles:

.. math:: F_i = \frac{m_i}{M} F

This weighting ensures that all the atoms experience the same acceleration, and hence move as one complete group.

There are three interaction available for dynamics run through *narupatools*:

Spring
------

The spring interaction is centered on a position :math:`p`, and has the following form:

.. math:: E = s \frac{k}{2} | p - r_{\text{COM}} |^2

.. math:: F = s k (p - r_{\text{COM}})

Here, :math:`k` is a scaling factor fixed as :math:`1 \text{kJ}/\text{mol}/\text{nm}^2`, and :math:`s` is the dimensionless scaling factor of the interaction.

Gaussian
--------

The gaussian interaction is centered on a position :math:`p`, and has the following form:

.. math:: E = s k \exp{\frac{| p - r_{\text{COM}} |^2}{2 \sigma^2}}

.. math:: F = s k (p - r_{\text{COM}}) \exp{\frac{| p - r_{\text{COM}} |^2}{2}}

Here, :math:`k` and :math:`\sigma` are scaling factors fixed at :math:`1 \text{kJ}/\text{mol}` and :math:`1 \text{nm}^2` respectively, and :math:`s` is the dimensionless scaling factor of the interaction.

Constant
--------

*narupatools* also adds a constant interaction, which applies a constant force :math:`F_c` to a set of particles:

.. math:: E = - F_c \cdot r_{\text{COM}}

.. math:: F = F_c

The dimensionless scaling factor is not applied to this kind of interaction.

Calculating Work
================

The work done by an interaction is given by the integral:

.. math:: W = \int_C F \cdot ds

where the integral C is along a given path.

Our trajectory is split into discrete steps, where we know both the force :math:`F` and displacement :math:`s` at different points in time.

The infinitesimal work here is given by:

.. math:: \delta W = F \cdot ds

For two adjacent timesteps with times :math:`t_0` and :math:`t_1`, we have:

.. math:: F = \frac{1}{2} (F(t_0) + F(t_1))

.. math:: ds = s(t_1) - s(t_0)

where we approximate the force as the average between the forces at the start and end of the step. Therefore:

.. math:: \delta W = \frac{1}{2} (F(t_0) + F(t_1)) \cdot (s(t_1) - s(t_0))

This gives the work for a single step on a single particle. The total work for an interaction is therefore:

.. math:: W = \frac{1}{2} \sum_i \sum_{n=0}^{N-1} (F_i(t_{n}) + F_i(t_{n+1})) \cdot (s_i(t_{n+1}) - s_i(t_{n}))

where the first sum is over all indices :math:`i` affected by the interaction and the second sum is over all timesteps