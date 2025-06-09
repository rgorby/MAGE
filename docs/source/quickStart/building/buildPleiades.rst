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

   module list
   No Modulefiles Currently Loaded.

Start by purging any currently-loaded modules, then loading the following
module set:

.. code-block:: bash

    module --force purge
    module load nas
    module load pkgsrc/2022Q1-rome
    module load comp-intel/2020.4.304
    module load mpi-hpe/mpt.2.23
    module load hdf5/1.8.18_mpt

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

    # Compile the kaiju software.
    make >& make.out

The build takes 20-30 minutes on ``pleiades``. When finished, your build
directory will contain a ``bin`` subdirectory which will contain the complete
set of ``kaiju`` executables.

.. note:: The compiled programs of interest in this case are
    ``voltron_mpi.x`` (for MAGE) and ``gamhelio_mpi.x`` (for the inner
    heliosphere). Documentation for the remaining tools is found
    :doc:`here </userGuide/kaijuTools/index>`. These tools are often run by
    using wrapper scripts found in the ``kaipy`` package.

.. note:: You can choose to compile only specific programs from the ``kaiju``
    package, which will speed up the compilation. For instance, to compile
    just the MPI components of ``MAGE``, you would use the command ``make
    gamera_mpi.x rcm.x remix.x voltron_mpi.x``.

