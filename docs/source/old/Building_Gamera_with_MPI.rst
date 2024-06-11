
Building Gamera with MPI
========================

The first step to building MPI versions of the gamera executables is making sure that you can build the non-MPI version of the gamera executables. Follow the [[Quick Start]], [[Gamerasphere]], and [[Gamera-RCM]] instructions to get to this point.

Once you are at the point that you can build the serial gamera executables, you need to modify your build environment to add MPI support. On Cheyenne this requires replacing the "hdf5/1.10.5" module with the "hdf5-mpi/1.10.5" module, and then also loading the "impi/2018.4.274" module. Your build environment should look like this:

.. code-block:: shell

   1) git/2.22.0
   2) intel/18.0.5
   3) hdf5-mpi/1.10.5
   4) impi/2018.4.274
   5) ncarenv/1.3
   6) ncarcompilers/0.5.0
   7) python/2.7.16
   8) cmake/3.14.4

On systems other than Cheyenne, your build environment may look different.

Additionally, you will need to enable the MPI option in cmake. You can do this by either adding this line to your kaiju/cmake/user.cmake file (creating this file if it does not exist):

.. code-block:: shell

   set(ENABLE_MPI ON CACHE BOOL "MPI" FORCE)

Or you can reconfigure your (or configure a new) cmake build folder with this command:

.. code-block:: shell

   cmake -DENABLE_MPI=ON ..

Either way, you will know that your build folder is configured correctly if a line like this is printed during the cmake configuration process:

.. code-block:: shell

   -- Found MPI: TRUE (found version "3.1") found components:  Fortran

Voltron MPI also uses MKL/pardiso which can be added to the build environment by loading the module 

.. code-block:: shell

   9) mkl/2018.0.5

and you can reconfigure cmake build folder with the commands

.. code-block:: shell

   ccmake .. (opens up cmake configuration)
   ENABLE_MKL ON
   [c], [g]

With that setup, you should be able to tell cmake to build the MPI executables:

.. code-block:: shell

   make gamera_mpi.x
   make voltron_mpi.x
