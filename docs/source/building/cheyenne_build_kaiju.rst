
Compiling the ``kaiju`` software on ``cheyenne``
========================================================

Introduction
------------

These instructions will walk you through the process of building and installing the ``kaiju`` software on the ``cheyenne.ucar.edu`` supercomputer.

These instructions assume that the user is using the ``bash`` shell, and that no modifications have been made to the user "dotfiles" (\ ``$HOME/.bashrc``\ , ``$HOME/.bash_profile``\ ). If you have customized either of these files for your account, please carefully save and inspect the output from each command in the build process to ensure that no unexpected problems have crept in. To facilitate this practice, all of the commands shown below will illustrate how to save command output, and how to measure how long each step takes. The latter is a useful bit of information which can help identify build problems early in the process, avoiding much wasted time and effort later.

WARNING
-------

The ``cheyenne`` and ``derecho`` systems share the same home directory for each user. However, these two systems contain different CPU architectures, and provide different module environments. Specifically, the default-provided module sets differ (significantly) between the two systems, as does the selection of available modules. Therefore, you *must* make sure that you keep software compiled for the two systems separate. Failure to do this will result in unpredictable (and worse, unreproducible) behavior on *both* systems.

In the instructions below, all code will be built and installed in a ``cheyenne``\ -specific subdirectory of the user home directory, i.e. ``$HOME/cheyenne``. This particular organization is not required - it is intended as an example of one possible way to segregate software that has been built for the two systems.

Building the serial version of the ``kaiju`` software on ``cheyenne``
-----------------------------------------------------------------------------

Step 1: Load the build modules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Begin by loading the modules needed to build the serial version of the ``kaiju`` software on ``cheyenne``\ :

.. code-block::

   #!shell
   module purge
   module load cmake/3.22.0
   module load git/2.33.1  # Needed for git-lfs
   module load ncarenv/1.3
   module load intel/2022.1
   module load geos/3.10.1  # Must come after intel/2022.1
   module load ncarcompilers/0.5.0  # Must come after intel/2022.1
   module load hdf5/1.12.2


After these commands have been run, verify your module list with ``module list``\ :

.. code-block::

   #!shell
   module list

   Currently Loaded Modules:
     1) cmake/3.22.0   2) git/2.33.1   3) ncarenv/1.3   4) intel/2022.1   5) geos/3.10.1   6) ncarcompilers/0.5.0   7) hdf5/1.12.2


NOTE: This list of modules will evolve. The current list is the recommended module set as of 2023-06-29.

Step 2: Create the build directory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Create a ``cheyenne``\ -specific build directory.

.. code-block::

   #!shell
   cd $HOME
   mkdir cheyenne
   cd cheyenne


Then make an additional subdirectory level for the branch of the code you are building (the ``development`` branch is used as an example). This arrangement is useful when you need to maintain simultaneous builds of different branches.

.. code-block::

   #!shell
   mkdir development
   cd development


Step 3: Clone the ``kaiju`` repository from BitBucket
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

NOTE: This step assumes you have been granted access to the ``kaiju`` repository on BitBucket, and that you have configured an SSH key pair for use with BitBucket. If you need help with these tasks, please contact a CGS team member for assistance.

Clone the ``kaiju`` repository (or "repo") from BitBucket:

.. code-block::

   #!shell
   git clone git@bitbucket.org:aplkaiju/kaiju.git


This process should take a minute or so. When complete, verify that the ``kaiju`` code exists in your directory (the actual directory contents may differ slightly from what is shown below):

.. code-block::

   #!shell
   ls kaiju
   analysis  cmake  CMakeLists.txt  examples  external  gitHookScripts  kaiju.sublime-project  kaipy  places  pytests  quickstart  README.md  scripts  setup.py  src  testingScripts  tests  xml


Now move down into the cloned repo, and switch to the branch of the code you wish to use. By default, the cloned repository provides the ``master`` branch, but we want the ``development`` branch:

.. code-block::

   #!shell
   cd kaiju
   git switch development


Step 4: Run ``cmake`` to create the ``Makefile`` needed to build the software
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Since the ``kaiju`` code can be built in serial and MPI forms, we first make a directory in which to build the serial version of the code (use whatever name you prefer, but ``build_serial`` is simple and unambiguous):

.. code-block::

   #!shell
   mkdir build_serial
   cd build_serial


Now run the ``cmake`` command. Save the ``cmake`` output, and use timestamps for each step. The options shown below direct the build process to use a recent version of the Intel Fortran compiler:

.. code-block::

   #!shell
   date; time FC=`which ifort` cmake -DALLOW_INVALID_COMPILERS=ON .. >& cmake.out; date


This command usually takes 10-20 seconds, depending on system activity. Examine the output file ``cmake.out`` for problems. It *should* look something like this:

.. code-block::

   -- The Fortran compiler identification is Intel 2021.5.0.20211109
   -- Detecting Fortran compiler ABI info
   -- Detecting Fortran compiler ABI info - done
   -- Check for working Fortran compiler: /glade/u/apps/ch/opt/ncarcompilers/0.5.0/intel/2022.1/ifort - skipped
   -- Found HDF5: Included by compiler wrappers  found components: Fortran 
   -- Found OpenMP_Fortran: -qopenmp (found version "5.0") 
   -- Found OpenMP: TRUE (found version "5.0") found components: Fortran 
   CMake Warning at cmake/compilers.cmake:61 (message):
     Setting default optimization to O2 to avoid certain Intel compiler bugs
   Call Stack (most recent call first):
     CMakeLists.txt:99 (include)


   -------------------------
   Configuration summary ...
   System: cheyenne4
   OS: Linux
   Processor: x86_64
   Compiler: Intel / 2021.5.0.20211109
             /glade/u/apps/ch/opt/ncarcompilers/0.5.0/intel/2022.1/ifort
   HDF5 Wrapper: 
   Version: 3f4e147c / development
   Build Type: Release
   Base Flags:  -fPIC -free -implicitnone -qopenmp
   Build Flags: -O2 -align array64byte -align rec32byte -no-prec-div -fast-transcendentals -ipo -march=corei7 -axCORE-AVX2
   -------------------------

   Adding CHIMP module ...
       EB IC file is /glade/u/home/ewinter/cheyenne/development/kaiju/src/chimp/ebICs/ebICstd.F90
       TP IC file is /glade/u/home/ewinter/cheyenne/development/kaiju/src/chimp/tpICs/tpICstd.F90
       Adding executable project.x
       Adding executable psd.x
       Adding executable push.x
       Adding executable slice.x
       Adding executable chop.x
       Adding executable trace.x
       Adding executable sctrack.x
       Adding executable calcdb.x
       Adding executable wpicheck.x
   Adding Gamera module ...
       Bricksize is 16
       IC file is /glade/u/home/ewinter/cheyenne/development/kaiju/src/gamera/ICs/null.F90
       Adding executable gamera.x
   Adding Gamera Helio module ...
       IC file is /glade/u/home/ewinter/cheyenne/development/kaiju/src/gamera/ICs/helio/wsa.F90
       Adding executable gamhelio.x
   Adding ReMIX module ...
       Adding executable remix.x
   Adding RCM module ...
       RCM Grid is of size 180 x 361 x 160
       Adding executable rcm.x
   Adding Voltron module ...
       IC file is /glade/u/home/ewinter/cheyenne/development/kaiju/src/voltron/ICs/earthcmi.F90
       Adding executable voltron.x
   -- Configuring done
   -- Generating done
   -- Build files have been written to: /glade/u/home/ewinter/cheyenne/development/kaiju/build_serial

Step 5: Compile the ``kaiju`` software
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now use ``make`` to build the ``kaiju`` software, time-stamping and saving the output:

.. code-block::

   #!shell
   date; time make >& make.out; date


This command should complete in about 20 minutes on ``cheyenne``. When the command is finished, check the output file ``make.out``. The file is long, but the last few lines should look something like this:

.. code-block::

   /glade/u/home/ewinter/cheyenne/cgs/kaiju/ewinter-python_environment_upgrade/kaiju/src/voltron/ICs/earthcmi.F90(169): remark #15009: uservolticinituser_mp_psphereic_ has been targeted for automatic cpu dispatch
   /glade/u/home/ewinter/cheyenne/cgs/kaiju/ewinter-python_environment_upgrade/kaiju/src/voltron/voltapp.F90(28): remark #15009: voltapp_mp_initvoltron_ has been targeted for automatic cpu dispatch
   /glade/u/home/ewinter/cheyenne/cgs/kaiju/ewinter-python_environment_upgrade/kaiju/src/base/kronos.F90(39): remark #15009: kronos_mp_initts_ has been targeted for automatic cpu dispatch
   [100%] Built target voltron.x
   [100%] Built target voltron

To verify that all of the ``kaiju`` programs have been built, examine the ``bin`` subdirectory of your ``build_serial`` directory (this list will evolve as more programs are added):

.. code-block::

   #!shell
   ls bin
   calcdb.x  chop.x  gamera.x  gamhelio.x  kaitoy.x  project.x  psd.x  push.x  rcm.x  remix2rcm.x  remix2remix.x  remix.x  sctrack.x  slice.x  trace.x  voltron.x  wpicheck.x



Building the MPI version of the ``kaiju`` software on ``cheyenne``
--------------------------------------------------------------------------

Step 1: Load the build modules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Begin by loading the modules needed to build the MPI version of the ``kaiju`` software on ``cheyenne``\ :

.. code-block::

   #!shell
   module purge
   module load cmake/3.22.0
   module load git/2.33.1  # Needed for git-lfs
   module load ncarenv/1.3
   module load intel/2022.1
   module load geos/3.10.1  # Must come after intel/2022.1
   module load ncarcompilers/0.5.0  # Must come after intel/2022.1
   module load mpt/2.25
   module load hdf5-mpi/1.12.2


After these commands have been run, verify your module list with ``module list``\ :

.. code-block::

   #!shell
   module list

   Currently Loaded Modules:
     1) cmake/3.22.0   2) git/2.33.1   3) ncarenv/1.3   4) intel/2022.1   5) geos/3.10.1   6) ncarcompilers/0.5.0   7) mpt/2.25   8) hdf5-mpi/1.12.2


NOTE: This list of modules will evolve. The current list is the recommended module set as of 2023-06-29.

Step 2: Create the build directory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Create a ``cheyenne``\ -specific build directory (skip if already done).

.. code-block::

   #!shell
   cd $HOME
   mkdir cheyenne
   cd cheyenne


Then make an additional subdirectory level for the branch of the code you are building (the ``development`` branch is used as an example). This arrangement is useful when you need to maintain simultaneous builds of different branches.

.. code-block::

   #!shell
   mkdir development
   cd development


Step 3: Clone the ``kaiju`` repository from BitBucket
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

NOTE: This step assumes you have been granted access to the ``kaiju`` repository on BitBucket, and that you have configured an SSH key pair for use with BitBucket. If you need help with these tasks, please contact a CGS team member for assistance.

Clone the ``kaiju`` repository (or "repo") from BitBucket (skip this part if you already cloned the repo):

.. code-block::

   #!shell
   git clone git@bitbucket.org:aplkaiju/kaiju.git


This process should take a minute or so. When complete, verify that the ``kaiju`` code exists in your directory (the actual directory contents may differ slightly from what is shown below):

.. code-block::

   #!shell
   ls kaiju
   analysis  cmake  CMakeLists.txt  examples  external  gitHookScripts  kaiju.sublime-project  kaipy  places  pytests  quickstart  README.md  scripts  setup.py  src  testingScripts  tests  xml


Now move down into the cloned repo, and switch to the branch of the code you wish to use. By default, the cloned repository provides the ``master`` branch, but we want the ``development`` branch:

.. code-block::

   #!shell
   cd kaiju
   git switch development


Step 4: Run ``cmake`` to create the ``Makefile`` needed to build the software
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Since the ``kaiju`` code can be built in serial and MPI forms, we first make a directory in which to build the MPI version of the code (use whatever name you prefer, but ``build_mpi`` is simple and unambiguous):

.. code-block::

   #!shell
   mkdir build_mpi
   cd build_mpi


Now run the ``cmake`` command. Save the ``cmake`` output, and use timestamps for each step. The options shown below direct the build process to use a recent version of the Intel Fortran compiler:

.. code-block::

   #!shell
   date; time FC=`which ifort` cmake -DENABLE_MPI=ON -DALLOW_INVALID_COMPILERS=ON .. >& cmake.out; date


This command usually takes 20-30 seconds, depending on system activity. Examine the output file ``cmake.out`` for problems. It *should* look something like this:

.. code-block::

   -- The Fortran compiler identification is Intel 2021.5.0.20211109
   -- Detecting Fortran compiler ABI info
   -- Detecting Fortran compiler ABI info - done
   -- Check for working Fortran compiler: /glade/u/apps/ch/opt/ncarcompilers/0.5.0/intel/2022.1/ifort - skipped
   Using default MPI Address Size of MPI_ADDRESS_KIND for basic MPI functions
   Using default MPI Address Size of MPI_ADDRESS_KIND for neighborhood MPI functions
   -- Found HDF5: /glade/u/apps/ch/opt/netcdf-mpi/4.9.0/mpt/2.25/intel/2022.1/lib/libhdf5_fortran.so;/glade/u/apps/ch/opt/netcdf-mpi/4.9.0/mpt/2.25/intel/2022.1/lib/libhdf5.so;/glade/u/apps/ch/opt/netcdf-mpi/4.9.0/mpt/2.25/intel/2022.1/lib/libsz.so;/glade/u/apps/ch/opt/netcdf-mpi/4.9.0/mpt/2.25/intel/2022.1/lib/libz.so;/usr/lib64/libdl.so;/usr/lib64/libm.so (found version "1.12.2") found components: Fortran 
   -- Found OpenMP_Fortran: -qopenmp (found version "5.0") 
   -- Found OpenMP: TRUE (found version "5.0") found components: Fortran 
   -- Found MPI
   CMake Warning at cmake/compilers.cmake:61 (message):
     Setting default optimization to O2 to avoid certain Intel compiler bugs
   Call Stack (most recent call first):
     CMakeLists.txt:99 (include)


   -------------------------
   Configuration summary ...
   System: cheyenne4
   OS: Linux
   Processor: x86_64
   Compiler: Intel / 2021.5.0.20211109
             mpif08
   HDF5 Wrapper: /glade/u/apps/ch/opt/netcdf-mpi/4.9.0/mpt/2.25/intel/2022.1/bin/h5pfc
   Version: 3f4e147c / development
   Build Type: Release
   Base Flags:  -fPIC -free -implicitnone -qopenmp -Wl,-Bstatic -Wl,-Bdynamic -Wl,--disable-new-dtags -Wl,-rpath,/glade/u/apps/ch/opt/geos/3.10.1/intel/2022.1/lib64 -Wl,-rpath,/glade/u/apps/ch/opt/netcdf-mpi/4.9.0/mpt/2.25/intel/2022.1/lib -L/glade/u/apps/ch/opt/geos/3.10.1/intel/2022.1/lib64 -lmpi
   Build Flags: -O2 -align array64byte -align rec32byte -no-prec-div -fast-transcendentals -ipo -march=corei7 -axCORE-AVX2
   -------------------------

   Adding CHIMP module ...
       EB IC file is /glade/u/home/ewinter/cheyenne/development/kaiju/src/chimp/ebICs/ebICstd.F90
       TP IC file is /glade/u/home/ewinter/cheyenne/development/kaiju/src/chimp/tpICs/tpICstd.F90
       Adding executable project.x
       Adding executable psd.x
       Adding executable push.x
       Adding executable slice.x
       Adding executable chop.x
       Adding executable trace.x
       Adding executable sctrack.x
       Adding executable calcdb.x
       Adding executable wpicheck.x
   Adding Gamera module ...
       Bricksize is 16
       IC file is /glade/u/home/ewinter/cheyenne/development/kaiju/src/gamera/ICs/null.F90
       Adding executable gamera.x
   Adding Gamera Helio module ...
       IC file is /glade/u/home/ewinter/cheyenne/development/kaiju/src/gamera/ICs/helio/wsa.F90
       Adding executable gamhelio.x
   Adding ReMIX module ...
       Adding executable remix.x
   Adding RCM module ...
       RCM Grid is of size 180 x 361 x 160
       Adding executable rcm.x
   Adding Voltron module ...
       IC file is /glade/u/home/ewinter/cheyenne/development/kaiju/src/voltron/ICs/earthcmi.F90
       Adding executable voltron.x
   Adding Base MPI module ...
       Adding executable kaitoy_mpi.x
   Adding Gamera MPI module ...
       Adding executable gamera_mpi.x
   Adding Voltron MPI module ...
       Adding executable voltron_mpi.x
   Adding Gamera Helio MPI module ...
       Adding executable gamhelio_mpi.x
   -- Configuring done
   -- Generating done
   -- Build files have been written to: /glade/u/home/ewinter/cheyenne/development/kaiju/build_mpi

Step 5: Compile the ``kaiju`` software
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now use ``make`` to build the ``kaiju`` software, time-stamping and saving the output:

.. code-block::

   #!shell
   date; time make >& make.out; date


This command should complete in about 25 minutes on ``cheyenne``. When the command is finished, check the output file ``make.out``. The file is long, but the last few lines should look something like this:

.. code-block::

   /glade/u/home/ewinter/cheyenne/development/kaiju/src/gamera/ICs/helio/wsa.F90(70): remark #15009: usergamic_mp_inituser_ has been targeted for automatic cpu dispatch
   /glade/u/home/ewinter/cheyenne/development/kaiju/src/gamera/ICs/helio/wsa.F90(222): remark #15009: usergamicinituser_mp_mapkinit_ has been targeted for automatic cpu dispatch
   /glade/u/home/ewinter/cheyenne/development/kaiju/src/gamera/ICs/helio/wsa.F90(201): remark #15009: usergamicinituser_mp_gasic_ has been targeted for automatic cpu dispatch
   [100%] Built target gamhelio_mpi.x
   [100%] Built target gamhelio_mpi

To verify that all of the ``kaiju`` programs have been built, examine the ``bin`` subdirectory of your ``build_mpi`` directory (this list will evolve as more programs are added):

.. code-block::

   #!shell
   ls bin
   calcdb.x  gamera_mpi.x  gamhelio_mpi.x  kaitoy_mpi.x  project.x  push.x  remix2rcm.x    remix.x    slice.x  voltron_mpi.x  wpicheck.x
   chop.x    gamera.x      gamhelio.x      kaitoy.x      psd.x      rcm.x   remix2remix.x  sctrack.x  trace.x  voltron.x


Using the ``kaiju`` software
--------------------------------

Once built, you must run the setup script before using the ``kaiju`` software:

.. code-block::

   #!shell
   source $KAIJU_HOME/scripts/setupEnvironment.sh


This script will set environment variables needed by the ``kaiju`` software, including the ``KAIJUHOME`` environment variable (not the ``KAIJU_HOME`` environment variable). However, the path to the compiled programs is not added - you will need to specify the complete path when using compiled programs. For example,  to run the serial version of ``gamera.x``\ :

.. code-block::

   #!shell
   $KAIJUHOME/build_serial/bin/gamera.x
