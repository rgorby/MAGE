Building the ``kaiju`` software on ``derecho``
==============================================

Introduction
------------

This page provides instructions for building the ``kaiju`` software on the
``derecho`` supercomputer.

Step 1: Preparing your software environment
-------------------------------------------

Like most HPC systems, ``derecho`` uses the ``module`` system to manage the
versions of software packages available to the user. A *module* is a
collection of programs and libraries for a specific task, and *loading*
the module adjusts the user environment (mostly by setting or updating
environment variables) to make that module available to the user. When you log
in to ``derecho``, the following modules are loaded by default:

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

Step 2: Build missing prerequisite libraries
--------------------------------------------

The `NASA CDF (Common Data Format) library <https://cdf.gsfc.nasa.gov/>`_ must
be built on ``derecho``. These instructions assume ``CDF_ROOT`` is the root of
the directory tree which will contain the CDF code and libraries. The build
follows the standard download/unpack/build/test/install cycle typical of most
open-source software.

.. important::

    The ``kaiju`` code has not been tested with versions of the CDF library
    more recent than the version indicated below.

.. code-block:: bash

   $ # Specify and create the root of the build tree.
   $ export CDF_ROOT=$HOME/local/cdf/3.9.0
   $ mkdir -p $CDF_ROOT/src
   $ cd $CDF_ROOT/src

   $ # Download the CDF source tarball.
   $ wget https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/cdf39_0/linux/cdf39_0-dist-all.tar.gz

   $ # Unpack the source code.
   $ tar xzvf cdf39_0-dist-all.tar.gz

   $ # Move into the code directory.
   $ cd cdf39_0-dist

   $ # Build the library using the default system compiler.
   $ make OS=linux ENV=gnu all

   $ # Test the library.
   $ make test

   $ # Install the library in a version-specific subdirectory.
   $ make INSTALLDIR=$CDF_ROOT install

   $ # Clean the build tree.
   $ make clean

   $ # Make the CDF library available for kaiju.
   $ source $CDF_ROOT/bin/definitions.B

Step 3: Create a ``python`` environment
---------------------------------------

Most of the ``kaiju`` software for pre-processing, post-processing, and
analysis is written in `Python <https://www.python.org/>`_. Python is available
in many forms (or 'distributions'), but we recommend use of the
`Miniconda distribution <https://docs.conda.io/en/latest/miniconda.html>`_ for
simplicity and compactness. These instructions install Miniconda into the home
directory. These instructions name the environment ``kaiju-3.8``, but you can
use any convenient name.

.. code-block:: bash

   $ # Download the Miniconda installer.
   $ cd $HOME
   $ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

   $ # Run the installer.
   $ bash ./Miniconda3-latest-Linux-x86_64.sh

   $ # Make sure the conda configuration is loaded.
   $ source $HOME/.bashrc

   $ # Create the environment for kaiju using python 3.8.
   $ conda create -n kaiju-3.8 python=3.8

   $ # Activate the new environment.
   $ conda activate kaiju-3.8

   $ # Install required packages.
   $ pip install alive_progress astropy cartopy cdasws configparser h5py matplotlib pandas progressbar pyhdf pyspedas pytest spacepy sunpy

Step 4: Build the ``kaiju`` software
------------------------------------

These instructions show how to build the MPI version of the ``kaiju``
software, since that is the version usually used on ``derecho``. These
instructions put the build in the subdirectory ``build_mpi`` under the
``kaiju`` source code directory, but you can place the build directory in any
convenient location.

.. code-block:: bash

    $ # Create the build directory.
    $ cd $HOME

    $ # Clone the ``kaiju`` repository and enter it.
    $ git clone git@bitbucket.org:aplkaiju/kaiju.git
    $ cd kaiju

    $ # Create the build directory and enter it.
    $ mkdir build_mpi
    $ cd build_mpi

    $ # Run cmake to create the Makefile.
    $ # The FC definition is required.
    $ FC=`which ifort` cmake ..

    $ # Compile the kaiju software.
    $ make

When finished, your build directory will contain a ``bin`` subdirectory which
should contain the following files (list current as of 4 March 2025):

.. code-block:: shell

    $ ls -w 80 bin
    calcdb.x        gamhelio.x    push.x         sctrack.x      voltron.x
    chop.x          kaitoy_mpi.x  rcm.x          slice.x        wpicheck.x
    gamera_mpi.x    kaitoy.x      remix2rcm.x    testNewRMS.x
    gamera.x        project.x     remix2remix.x  trace.x
    gamhelio_mpi.x  psd.x         remix.x        voltron_mpi.x