Systems
=======

The particles that we our running dynamics on will be referred to as a *system*. This is a brief exposition on some of the concepts involved in this.

A typical molecular dynamics system consists of one or more particles. Commonly, these particles may represent atoms, but there is no requirement that they do. At the very least, a particle has a **position** in three dimensional space and a **mass**. If dynamics are being run, then the particles will also have a **velocity** which dictates how quickly they are moving through space and at what rate. In turn, the velocity is changed by any **acceleration** that is applied each atom. For a classical system, this acceleration occurs as a result a **force** applied to each atom, due to interactions with each other. The two are related by **Newton's Second Law**:

.. math:`F_i = m_i a_i`

From this, it's clear to see that the mass in the context of this equation represents the particle's resistance to change. Namely, if a particle has a higher mass it requires more force to cause the same acceleration.