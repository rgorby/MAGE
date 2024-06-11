
This is a simple set of instructions to run an MHD test case w/ Gamera through the Kaiju repo.

Installation
^^^^^^^^^^^^

Use git clone to download the repository to where you want it.  

The aplkaiju bitbucket site has the current version of the coupled code base.

Assuming $KAIJUDIR is base installation directory, add $KAIJUDIR and $KAIJUDIR/scripts to $PYTHONPATH, and just $KAIJUDIR/scripts to $PATH.

Loading Modules
^^^^^^^^^^^^^^^

Load the following modules in the following order on Pleiades:

.. code-block::

   module load pkgsrc
   module load comp-intel/2020.4.304
   module load mpi-hpe/mpt.2.23
   module load hdf5/1.8.18_mpt
   module load python3/3.7.0

Building
^^^^^^^^

With the modules loaded, you should be able to build the MPI gamera and voltron applications. Navigate into the the root repository folder ($KAIJUDIR above), create a new folder for the build process, and enter this folder. This folder can be named anything, I will call it "build" for this example.

.. code-block::

   mkdir build
   cd build

Once in this folder, call cmake on the root of the repository to perform the cmake process. Be sure to enable MPI with the ENABLE_MPI flag. You can then build the gamera_mpi and voltron_mpi executables.

.. code-block::

   cmake -DENABLE_MPI=ON ..
   make gamera_mpi.x
   make voltron_mpi.x

Using MPT
^^^^^^^^^

Because Pleiades uses the MPT mpi library, the submission scripts will require modifications compared to the ones typically used on Cheyenne. These changes are described in a dedicated wiki page [[Running with MPT]].

Notes
=====

The build folder does not have to be located inside the repository folder, but it is common to do so. If you locate it somewhere else, just be sure to adjust the cmake command so that it points to the root of the repository.

There are other MPI modules on Pleiades, such as Intel-MPI modules located at /nasa/modulefiles/testing/mpi-intel/ . These are not tested, and use-at-your-own-risk. Using them in place of the MPT MPI module will require replacing other modules, as well.
