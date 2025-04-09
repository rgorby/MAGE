Building the ``kaiju`` software
===============================

Introduction
------------

This section describes how to build the ``kaiju`` software on two different
supercomputers - ``derecho`` and ``pleiades``. If you are trying to build the
``kaiju`` software on a different system, use these instructions as a starting
point.

Before you begin
----------------

The ``kaiju`` software is typically built with the Intel Fortran compiler
(although the software can also be built with any Fortran compiler,
e.g., GNU).
Building the software also requires the ``cmake`` build tool, the ``HDF5``
library, and an ``MPI`` library. Instructions for loading these packages are
provided in the ``module`` commands for each HPC system.

Getting the source code
-----------------------

The ``kaiju`` source code can be obtained by cloning the ``kaiju`` repository
on BitBucket:

.. code-block:: bash

    $ git clone https://bitbucket.org/aplkaiju/kaiju.git

.. important::

    The ``kaiju`` repository on BitBucket uses ``git-lfs`` to support the
    use of large binary files. You *must* make sure ``git-lfs`` is
    available in your ``git`` installation to ensure a complete clone
    of the ``kaiju`` repository.

Table of Contents
-----------------

.. toctree::
    :maxdepth: 1

    buildDerecho
    buildPleiades
