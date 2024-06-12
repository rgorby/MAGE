
This is a simple set of instructions to run an MHD test case w/ Gamera through the Kaiju repo.

Installation
============

Use git clone to download the repository to where you want it.

The aplkaiju bitbucket site has the current version of the coupled code base.

Assuming $KAIJUDIR is base installation directory, add $KAIJUDIR to $PYTHONPATH, and $KAIJUDIR/scripts to $PATH.

Loading Modules
===============

Load the following modules on Frontera in addition to those that are pre-loaded:

.. code-block:: python

   module load intel/19.1.3
   module load phdf5

this will give a module setup that looks like below

.. code-block:: python

     1) git/2.24.1      4) pmix/3.1.4      7) TACC          10) python3/3.7.0
     2) autotools/1.2   5) hwloc/1.11.12   8) intel/19.1.3  11) phdf5/1.10.4
     3) cmake/3.20.3    6) xalt/2.10.31    9) impi/19.0.9

As of the writing of this the intel/19.1.3 module was not publicly available on Frontera so you will have to download it. 

Building
========

In the folder $KAIJUDIR/cmake create a file "user.cmake" to set some compiler and linker flags as below:

.. code-block:: python

   string(APPEND CMAKE_Fortran_FLAGS " -I$ENV{TACC_HDF5_INC}")
   string(APPEND CMAKE_Fortran_FLAGS " -L$ENV{TACC_HDF5_LIB} -lhdf5_fortran -lhdf5hl_fortran")

With modules loaded and flags set, you should be able to build the MPI gamera and voltron applications. Navigate into the the root repository folder ($KAIJUDIR above), create a new folder for the build process, and enter this folder. This folder can be named anything, I will call it "build" for this example.

.. code-block:: python

   mkdir build
   cd build

Once in this folder, set the $FC environment variable and call cmake on the root of the repository to perform the cmake process. Be sure to enable MPI with the ENABLE_MPI flag. You can then build the gamera_mpi and voltron_mpi executables.

.. code-block:: python

   export FC="h5pfc"
   cmake -DENABLE_MPI=ON ..
   make gamera_mpi.x
   make voltron_mpi.x
