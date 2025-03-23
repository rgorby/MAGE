Building the ``kaiju`` software on ``pleiades``
===============================================

Introduction
------------

This page provides instructions for building the ``kaiju`` software on the
``pleiades`` supercomputer. These instructions assume that you have cloned the
``kaiju`` repository.

Prepare your software environment
---------------------------------

Like most HPC systems, ``pleiades`` uses the ``module`` system to manage the
versions of software packages available to the user. When you log in to
``pleiades``, no modules are loaded by default:

.. code-block:: bash

   $ module list
   No Modulefiles Currently Loaded.

Start by purging any currently-loaded modules, then loading the following
modules:

.. code-block:: bash

    $ module --force purge
    $ module load nas
    $ module load pkgsrc/2022Q1-rome
    $ module load comp-intel/2020.4.304
    $ module load mpi-hpe/mpt.2.23
    $ module load hdf5/1.8.18_mpt

.. important::

    You must use these exact versions of the modules to ensure the software
    compiles properly. If you use different versions of any of these modules,
    a successful build cannot be guaranteed. This module list is current as of
    **4 March 2025**, and is subject to change as the compute environment
    changes.

Build the ``kaiju`` software
----------------------------

These instructions show how to build the MPI version of the ``kaiju``
software, since that is the version usually used on ``pleiades``. These
instructions put the build in the subdirectory ``build_mpi`` under the
``kaiju`` source code directory, but you can place the build directory in any
convenient location.

.. code-block:: bash

    $ # Move to your kaiju clone.
    cd /path/to/kaiju

    $ # Create the build directory and enter it.
    $ mkdir build_mpi
    $ cd build_mpi

    $ # Run cmake to create the Makefile, saving output.
    $ # The FC definition is required.
    $ FC=`which ifort` cmake -DENABLE_MPI=ON .. >& cmake.out

    $ # Compile the kaiju software.
    $ make >& make.out

The build takes 20-30 minutes on ``pleiades``. When finished, your build
directory will contain a ``bin`` subdirectory which should contain the
following files (list current as of 4 March 2025):

.. code-block:: shell

    $ ls -w 80 bin
    calcdb.x        gamhelio.x    push.x         sctrack.x      voltron.x
    chop.x          kaitoy_mpi.x  rcm.x          slice.x        wpicheck.x
    gamera_mpi.x    kaitoy.x      remix2rcm.x    testNewRMS.x
    gamera.x        project.x     remix2remix.x  trace.x
    gamhelio_mpi.x  psd.x         remix.x        voltron_mpi.x