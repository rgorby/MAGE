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

   $ module list

   Currently Loaded Modules:
     1) ncarenv/23.09  (S)   4) ncarcompilers/1.0.0   7) netcdf/4.9.2
     2) craype/2.7.23        5) cray-mpich/8.1.27
     3) intel/2023.2.1       6) hdf5/1.12.2

     Where:
      S:  Module is Sticky, requires --force to unload or purge

This set of modules **cannot** be used to build the ``kaiju`` code.

Start by purging any currently-loaded modules, then loading the following
modules:

.. code-block:: bash

    $ module --force purge
    $ module load ncarenv/23.06
    $ module load cmake/3.26.3
    $ module load craype/2.7.20
    $ module load intel/2023.0.0
    $ module load ncarcompilers/1.0.0
    $ module load cray-mpich/8.1.25
    $ module load hdf5-mpi/1.12.2

.. important::

    You must use these exact versions of the modules to ensure the software
    compiles properly. If you use different versions of any of these modules,
    a successful build cannot be guaranteed. This module list is current as of
    **4 March 2025**, and is subject to change as the compute environment
    changes.

Build the ``kaiju`` software
----------------------------

These instructions show how to build the MPI version of the ``kaiju``
software, since that is the version usually used on ``derecho``. These
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

The build takes about 15 minutes on ``derecho``. When finished, your build
directory will contain a ``bin`` subdirectory which should contain the
following files (list current as of 4 March 2025):

.. attention:: 

               please add instructions for non-MPI builds (voltron.x
               and gamera.x)

.. code-block:: shell

    $ ls -w 80 bin
    calcdb.x        gamhelio.x    push.x         sctrack.x      voltron.x
    chop.x          kaitoy_mpi.x  rcm.x          slice.x        wpicheck.x
    gamera_mpi.x    kaitoy.x      remix2rcm.x    testNewRMS.x
    gamera.x        project.x     remix2remix.x  trace.x
    gamhelio_mpi.x  psd.x         remix.x        voltron_mpi.x


.. attention:: 
               are we describing these tools in this documentation
               (e.g., calcdb.x and sctrack.x? if so, please link below)
    
.. note:: The compiled programs of interest in this case are ``voltron_mpi.x``
    (used for ``MAGE``) and ``gamhelio_mpi.x``
    (used for the inner heliosphere). The tool ``calcdb.x`` computes ground
    delta-B values from the results of a MAGE model, and ``sctrack.x``
    interpolates model results to spacecraft trajectories. These tools, and
    most of the remaining compiled programs, are typically run by using
    wrapper scripts found in the ``kaipy`` package. Documentation for
    ``kaipy`` is available
    `here <https://kaipy-docs.readthedocs.io/en/latest/index.html>`_.

.. note:: Note also that you can choose to compile only the target
          that you are interested in, which will speed up the
          compilation. For instance, to compile ``MAGE``, you would do

.. attention::
   Eric, please, check and complete below 
          
.. code-block::

   $ make voltron_mpi.x
