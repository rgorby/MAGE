Building the ``kaiju`` software on ``derecho``
==============================================


Introduction
------------

This page provides instructions for building the ``kaiju`` software on the
``derecho`` supercomputer. These instructions assume that you have cloned the
``kaiju`` repository.


Prepare your software environment
---------------------------------

Like most HPC systems, ``derecho`` uses the ``module`` system to manage the
versions of software packages available to the user. When you log in to
``derecho``, the following modules are loaded by default:

.. code-block:: bash

   module list

   Currently Loaded Modules:
     1) ncarenv/23.09  (S)   4) ncarcompilers/1.0.0   7) netcdf/4.9.2
     2) craype/2.7.23        5) cray-mpich/8.1.27
     3) intel/2023.2.1       6) hdf5/1.12.2

     Where:
      S:  Module is Sticky, requires --force to unload or purge

This set of modules **cannot** be used to build the ``kaiju`` code.

Start by purging any currently-loaded modules, then loading the following
module set:

.. code-block:: bash

    module --force purge
    module load ncarenv/23.06
    module load cmake/3.26.3
    module load craype/2.7.20
    module load intel/2023.0.0
    module load ncarcompilers/1.0.0
    module load cray-mpich/8.1.25
    module load hdf5-mpi/1.12.2

.. important::

    You must use these exact versions of the modules to ensure the software
    compiles properly. If you use different versions of any of these modules,
    a successful build cannot be guaranteed. This module list is current as of
    **11 April 2025**, and is subject to change as the compute environment
    changes.


Build the ``kaiju`` software
----------------------------

These instructions show how to build the MPI version of the ``kaiju``
software. The MPI version is built in the subdirectory ``build_mpi``
under the ``kaiju`` source code directory. In practice, you can place the
build directory in any convenient location.

.. code-block:: bash

    # Move to your kaiju clone.
    cd /path/to/kaiju

    # Create the build directory and enter it.
    mkdir build_mpi
    cd build_mpi

    # Run cmake to create the Makefile, saving output.
    # NOTE: The FC definition is *required* for proper cmake operation.
    FC=`which ifort` cmake -DENABLE_MPI=ON .. >& cmake.out

    $ # Compile the kaiju software.
    $ make >& make.out

The MPI build takes about 15 minutes on ``derecho``. When finished, your build
directory will contain a ``bin`` subdirectory which will contain the complete
set of ``kaiju`` executables.

.. code-block:: shell

    ls -w 80 bin
    calcdb.x        gamhelio.x    push.x         sctrack.x      voltron.x
    chop.x          kaitoy_mpi.x  rcm.x          slice.x        wpicheck.x
    gamera_mpi.x    kaitoy.x      remix2rcm.x    testNewRMS.x
    gamera.x        project.x     remix2remix.x  trace.x
    gamhelio_mpi.x  psd.x         remix.x        voltron_mpi.x

.. note:: The compiled programs of interest in this case are
    ``voltron_mpi.x``, ``gamera_mpi.x``, ``rcm.x`` and ``remix.x`` (all are
    part of the MAGE model),  and ``gamhelio_mpi.x`` (used for the inner
    heliosphere). Documentation for the most useful of the remaining tools is
    found :doc:`here </userGuide/kaijuTools/index>`. These tools are often
    run by using wrapper scripts found in the ``kaipy`` package. Documentation
    for ``kaipy`` is available
    `here <https://kaipy-docs.readthedocs.io/en/latest/index.html>`_.

.. note:: You can choose to compile only specific programs from the ``kaiju``
    package, which will speed up the compilation. For instance, to compile
    just the MPI components of ``MAGE``, you would use the command
    ``make gamera_mpi.x rcm.x remix.x voltron_mpi.x``.
