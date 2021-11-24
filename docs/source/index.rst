narupatools
===========

.. image:: https://img.shields.io/badge/Source-GitHub-181717?style=flat-square&logo=github
    :target: https://github.com/alexjbinnie/narupatools
.. image:: https://img.shields.io/badge/Download-Anaconda-44A833?style=flat-square&logo=anaconda
    :target: https://anaconda.org/alexjbinnie/narupatools
.. image:: https://img.shields.io/badge/Documentation-Read_the_Docs-8CA1AF?style=flat-square&logo=read-the-docs
    :target: https://narupatools.readthedocs.io/en/latest/
.. image:: https://img.shields.io/badge/License-GNU_General_Public_License_V3-orange?style=flat-square
    :target: https://www.gnu.org/licenses/gpl-3.0.en.html

*narupatools* is a python package containing utilities for running and interacting with molecular dynamics simulations being run within
the `Narupa framework <https://gitlab.com/intangiblerealities/narupa-protocol>`_.

This documentation contains an overview of some of the concepts involved in *narupatools*, tutorials describing how to use the various parts, and documentation for all public methods and classes.

Installation
------------

.. image:: https://anaconda.org/alexjbinnie/narupatools/badges/version.svg
    :target: https://anaconda.org/alexjbinnie/narupatools
.. image:: https://anaconda.org/alexjbinnie/narupatools/badges/latest_release_date.svg
    :target: https://anaconda.org/alexjbinnie/narupatools

You can install *narupatools* using conda:

.. code-block:: console

   conda install -c conda-forge -c irl -c alexjbinnie narupatools

*narupatools* has Narupa and its dependencies as dependencies itself, so if setting up a new conda environment you can simply use:

.. code-block:: console

   conda create -n narupatools
   conda activate narupatools
   conda install-c conda-forge -c irl -c alexjbinnie narupatools

Disclaimer
----------

*narupatools* is a separate project which is not a part of Narupa, and does not have any affiliation with either Narupa or the Intangible Realities Laboratory. The author would like to make it clear that issues arising due to this library should not be reported on the Narupa repositories, and instead reported to the narupatools `bug tracker <https://gitlab.com/alexjbinnie/narupatools/-/issues>`_.

License
-------

This project is released under the GNU General Public License version 3. All code written for narupatools is (c) Alex Jamieson-Binnie, All Rights Reserved.

Some code is originally part of the narupa python libraries, and is modified under the terms of the GPL. All narupa code is (c) Intangible Realities Lab, University Of Bristol. All rights reserved.

.. toctree::
   :maxdepth: 2

   FAQ <faq/faq>
   Concepts <concepts/concepts>
   Tutorials <tutorials/tutorials>

   API Reference <api/narupatools>


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
