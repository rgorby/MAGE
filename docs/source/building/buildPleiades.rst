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
    
    module use -a /swbuild/analytix/tools/modulefile
    module load nas
    module load pkgsrc/2022Q1-rome
    module load comp-intel/2020.4.304
    module load mpi-hpe/mpt.2.23
    module load hdf5/1.8.18_mpt
    module load miniconda3/v4


.. note::

    For MAGE runs coupled with TIEGCM (known as "GTR"), use the following modules:

    .. warning::

        The GTR currently required custom built NetCDF and ESMF modules on ``pleiades``. If you need to
        run GTR, you will need access to ``/home7/nrao3/local3`` and ``/nobackup/nrao3/tiegcm/tiegcm3.0/data``,
        please reach out to ``nikhilr@ucar.edu`` with the following:
        
            - Your Pleiades username
            - Your Name
            - Your Institution

    .. code-block:: bash

        module --force purge
        
        module use -a /nasa/modulefiles/testing
        module use -a /swbuild/analytix/tools/modulefile
        module load nas
        module load comp-intel/2020.4.304
        module load mpi-hpe/mpt.2.30
        module load szip/2.1.1
        module load hdf5/1.12.3_mpt
        module load miniconda3/v4

        export FC=mpif90
        export CC=mpicc
        export CXX=mpicxx

        export PREFIX=/home7/nrao3/local3
        export LIBRARY_PATH=${LIBRARY_PATH}:$PREFIX/lib
        export LD_LIBRARY_PATH=$LIBRARY_PATH
        export CPATH=$PREFIX/include
        export PATH=${PATH}:$PREFIX/bin
    
    .. warning::
        
        Build the GTR version of the code in the ``build_gtr`` subdirectory
        under the ``kaiju`` source code directory, not the ``build_mpi``
        subdirectory as the module set is different.


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

    # If you are building the GTR version of the code, create the build_gtr directory instead:
    # mkdir build_gtr
    # cd build_gtr


    # Run cmake to create the Makefile, saving output.
    # NOTE: The FC definition is *required* for proper cmake operation.
    FC=`which ifort` cmake -DENABLE_MPI=ON .. >& cmake.out

    # You can pick one compile target below or compile all of them, if you'd like

    # Compile the MAGE model for geospace simulations
    make -j4 voltron_mpi.x >& make-voltron.out

    # Compile the GAMERA-helio model for inner heliosphere simulations
    make -j4 gamhelio_mpi.x >& make-gamhelio.out

    # Compile analysis tools
    make -j4 calcdb.x chop.x sctrack.x slice.x >& make-analysis.out

    
When finished, your build directory will contain a ``bin``
subdirectory which will contain the compiled ``kaiju`` executables.

.. note:: Documentation on the analysis tools is found
    :doc:`here </tools/index>`.


