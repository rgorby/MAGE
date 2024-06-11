
Building the ``kaiju`` software on ``derecho``
======================================================

Introduction
------------

These instructions will walk you through the process of building and installing the ``kaiju`` software on the ``derecho`` supercomputer.

These instructions assume that the user is using the ``bash`` shell, and that no modifications have been made to the user "dotfiles" (\ ``$HOME/.bashrc``\ , ``$HOME/.bash_profile``\ ). If you have customized either of these files for your account, please carefully save and inspect the output from each command in the build process to ensure that no unexpected problems have crept in. To facilitate this practice, all of the commands shown below will illustrate how to save command output, and how to measure how long each step takes. The latter is a useful bit of information which can help identify build problems early in the process, avoiding much wasted time and effort later.

Like most HPC systems, ``derecho`` uses the ``module`` system to manage the versions of software packages available to the user. A "module" is a collection of programs and libraries for a specific task, and "loading" the module adjusts the user environment (mostly by setting or updating environment variables) to make that module available to the user. When you log in to ``derecho``\ , the following modules are loaded by default:

.. code-block::

   #!shell
   module list

   Currently Loaded Modules:
     1) ncarenv/23.06 (S)   2) craype/2.7.20   3) intel/2023.0.0   4) ncarcompilers/1.0.0   5) cray-mpich/8.1.25   6) hdf5/1.12.2   7) netcdf/4.9.2

      Where:
       S:  Module is Sticky, requires --force to unload or purge


This default set of modules needs a few changes in order to build the ``kaiju`` software. The correct module sets for the serial and MPI versions of ``kaiju`` will be listed in the instructions below.

NOTE: In the commands shown below, the shell prompt (the strings starting with ``ewinter@derecho4:``\ ) are provided to distinguish commands from their outputs. In your session, your prompt may be different.

WARNING
-------

The ``casper``\ , ``cheyenne`` and ``derecho`` systems share the same home directory for each user. However, these three systems contain different CPU architectures, and provide different module environments. Specifically, the default module sets differ (significantly) between the three systems, as does the selection of available modules. Therefore, you *must* make sure that you keep software compiled for the three systems separate. Failure to do this will result in unpredictable (and worse, unreproducible) behavior on *all three* systems.

In the instructions below, all code will be built and installed in a ``derecho``\ -specific subdirectory of the user home directory, i.e. ``$HOME/derecho``. This particular organization is not required - it is intended as an example of one possible way to segregate software that has been built for the three systems.

Building the serial version of the ``kaiju`` software on ``derecho``
----------------------------------------------------------------------------

Step 1: Load the build modules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Begin by loading the modules needed to build the serial version of the ``kaiju`` software on ``derecho``\ :

.. code-block::

   #!shell
   module purge
   module load ncarenv/23.06
   module load craype/2.7.20
   module load intel/2023.0.0
   module load ncarcompilers/1.0.0  # Must come after intel/2023.0.0
   module load hdf5/1.12.2
   module load cmake/3.26.3
   module load geos/3.9.1  # Must come after intel/2023.0.0
   module list

After these commands have been run, verify your module list with ``module list``\ :

.. code-block::

   #!shell
   module list

   Currently Loaded Modules:
   1) ncarenv/23.06 (S)   3) intel/2023.0.0   5) ncarcompilers/1.0.0   7) cmake/3.26.3
   2) craype/2.7.20       4) geos/3.9.1       6) hdf5/1.12.2


NOTE: This list of modules will evolve. The current list is the recommended module set as of 2023-06-29.

Step 2: Create the build directory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Create a ``derecho``\ -specific build directory.

.. code-block::

   #!shell
   cd $HOME
   mkdir derecho
   cd derecho


Then make an additional subdirectory level for the branch of the code you are building (the ``development`` branch is used as an example). This arrangement is useful when you need to maintain simultaneous builds of different branches.

.. code-block::

   #!shell
   mkdir development
   development


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
   analysis        examples        kaiju.sublime-project  pytests     scripts   testingScripts
   cmake           external        kaipy                  quickstart  setup.py  tests
   CMakeLists.txt  gitHookScripts  places                 README.md   src       xml


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


This command usually takes about 5 seconds, depending on system activity. Examine the output file ``cmake.out`` for problems. It *should* look something like this:

.. code-block::

   -- The Fortran compiler identification is Intel 2021.8.0.20221119
   -- Cray Programming Environment 2.7.20 Fortran
   -- Detecting Fortran compiler ABI info
   -- Detecting Fortran compiler ABI info - done
   -- Check for working Fortran compiler: /glade/u/apps/derecho/23.06/spack/opt/spack/ncarcompilers/1.0.0/oneapi/2023.0.0/ec7b/bin/ifort - skipped
   -- Found HDF5: hdf5_fortran-shared (found version "1.12.2") found components: Fortran 
   -- Found OpenMP_Fortran: -qopenmp (found version "5.0") 
   -- Found OpenMP: TRUE (found version "5.0") found components: Fortran 
   CMake Warning at cmake/compilers.cmake:61 (message):
     Setting default optimization to O2 to avoid certain Intel compiler bugs
   Call Stack (most recent call first):
     CMakeLists.txt:99 (include)


   -------------------------
   Configuration summary ...
   System: derecho4
   OS: Linux
   Processor: x86_64
   Compiler: Intel / 2021.8.0.20221119
             /glade/u/apps/derecho/23.06/spack/opt/spack/ncarcompilers/1.0.0/oneapi/2023.0.0/ec7b/bin/ifort
   HDF5 Wrapper: 
   Version: 3f4e147c / development
   Build Type: Release
   Base Flags:  -fPIC -free -implicitnone -qopenmp
   Build Flags: -O2 -align array64byte -align rec32byte -no-prec-div -fast-transcendentals -ipo
   -------------------------

   Adding CHIMP module ...
       EB IC file is /glade/u/home/ewinter/derecho/development/kaiju/src/chimp/ebICs/ebICstd.F90
       TP IC file is /glade/u/home/ewinter/derecho/development/kaiju/src/chimp/tpICs/tpICstd.F90
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
       IC file is /glade/u/home/ewinter/derecho/development/kaiju/src/gamera/ICs/null.F90
       Adding executable gamera.x
   Adding Gamera Helio module ...
       IC file is /glade/u/home/ewinter/derecho/development/kaiju/src/gamera/ICs/helio/wsa.F90
       Adding executable gamhelio.x
   Adding ReMIX module ...
       Adding executable remix.x
   Adding RCM module ...
       RCM Grid is of size 180 x 361 x 160
       Adding executable rcm.x
   Adding Voltron module ...
       IC file is /glade/u/home/ewinter/derecho/development/kaiju/src/voltron/ICs/earthcmi.F90
       Adding executable voltron.x
   -- Configuring done (3.9s)
   -- Generating done (0.9s)
   -- Build files have been written to: /glade/u/home/ewinter/derecho/development/kaiju/build_serial

Step 5: Compile the ``kaiju`` software
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now use ``make`` to build the ``kaiju`` software, time-stamping and saving the output:

.. code-block::

   #!shell
   date; time make >& make.out; date


This command should complete in about 6-7 minutes on ``derecho`` (yes, ``derecho`` is much faster than ``cheyenne``\ ). When the command is finished, check the output file ``make.out``. The file is long, but the last few lines should look something like this:

.. code-block::

   [ 98%] Building Fortran object src/voltron/CMakeFiles/voltlib.dir/voltio.F90.o
   [ 99%] Building Fortran object src/voltron/CMakeFiles/voltlib.dir/voltapp.F90.o
   [ 99%] Linking Fortran static library libvoltlib.a
   [ 99%] Built target voltlib
   [100%] Building Fortran object CMakeFiles/voltron.x.dir/src/drivers/voltronx.F90.o
   [100%] Linking Fortran executable bin/voltron.x
   [100%] Built target voltron.x
   [100%] Built target voltron

To verify that all of the ``kaiju`` programs have been built, examine the ``bin`` subdirectory of your ``build_serial`` directory (this list will evolve as more programs are added):

.. code-block::

   #!shell
   ls bin
   calcdb.x  gamera.x    kaitoy.x   psd.x   rcm.x        remix2remix.x  sctrack.x  trace.x    wpicheck.x
   chop.x    gamhelio.x  project.x  push.x  remix2rcm.x  remix.x        slice.x    voltron.x



Building the MPI version of the ``kaiju`` software on ``derecho.hpc.ucar.edu``
--------------------------------------------------------------------------------------

Step 1: Load the build modules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Begin by loading the modules needed to build the serial version of the ``kaiju`` software on ``derecho``\ :

.. code-block::

   #!shell
   module purge
   module load cmake/3.26.3
   module load craype/2.7.20
   module load intel/2023.0.0
   module load geos/3.9.1  # Must come after intel/2023.0.0
   module load ncarcompilers/1.0.0  # Must come after intel/2023.0.0
   module load cray-mpich/8.1.25
   module load hdf5-mpi/1.12.2


After these commands have been run, verify your module list with ``module list``\ :

.. code-block::

   #!shell
   module list

   Currently Loaded Modules:
     1) ncarenv/23.06 (S)   3) craype/2.7.20    5) geos/3.9.1            7) cray-mpich/8.1.25
     2) cmake/3.26.3        4) intel/2023.0.0   6) ncarcompilers/1.0.0   8) hdf5-mpi/1.12.2


NOTE: This list of modules will evolve. The current list is the recommended module set as of 2023-06-29.

Step 2: Create the build directory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Create a ``derecho``\ -specific build directory (skip if already done).

.. code-block::

   #!shell
   cd $HOME
   mkdir derecho
   cd derecho


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
   analysis      CMakeLists.txt  gitHookScripts         places      README.md  src             xml
   build_serial  examples        kaiju.sublime-project  pytests     scripts    testingScripts
   cmake         external        kaipy                  quickstart  setup.py   tests


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


This command usually takes 12-15 seconds, depending on system activity. Examine the output file ``cmake.out`` for problems. It *should* look something like this:

.. code-block::

   -- The Fortran compiler identification is Intel 2021.8.0.20221119
   -- Cray Programming Environment 2.7.20 Fortran
   -- Detecting Fortran compiler ABI info
   -- Detecting Fortran compiler ABI info - done
   -- Check for working Fortran compiler: /glade/u/apps/derecho/23.06/spack/opt/spack/ncarcompilers/1.0.0/oneapi/2023.0.0/ec7b/bin/ifort - skipped
   Using default MPI Address Size of MPI_ADDRESS_KIND for basic MPI functions
   Using default MPI Address Size of MPI_ADDRESS_KIND for neighborhood MPI functions
   -- Found HDF5: hdf5_fortran-shared (found version "1.12.2") found components: Fortran 
   -- Found OpenMP_Fortran: -qopenmp (found version "5.0") 
   -- Found OpenMP: TRUE (found version "5.0") found components: Fortran 
   -- Found MPI_Fortran: /glade/u/apps/derecho/23.06/spack/opt/spack/ncarcompilers/1.0.0/oneapi/2023.0.0/ec7b/bin/ifort (found version "3.1") 
   -- Found MPI: TRUE (found version "3.1") found components: Fortran 
   -- Found MPI
   CMake Warning at cmake/compilers.cmake:61 (message):
     Setting default optimization to O2 to avoid certain Intel compiler bugs
   Call Stack (most recent call first):
     CMakeLists.txt:99 (include)


   -------------------------
   Configuration summary ...
   System: derecho4
   OS: Linux
   Processor: x86_64
   Compiler: Intel / 2021.8.0.20221119
             /glade/u/apps/derecho/23.06/spack/opt/spack/ncarcompilers/1.0.0/oneapi/2023.0.0/ec7b/bin/ifort
   HDF5 Wrapper: 
   Version: 3f4e147c / development
   Build Type: Release
   Base Flags:  -fPIC -free -implicitnone -qopenmp  -lmpi
   Build Flags: -O2 -align array64byte -align rec32byte -no-prec-div -fast-transcendentals -ipo
   -------------------------

   Adding CHIMP module ...
       EB IC file is /glade/u/home/ewinter/derecho/development/kaiju/src/chimp/ebICs/ebICstd.F90
       TP IC file is /glade/u/home/ewinter/derecho/development/kaiju/src/chimp/tpICs/tpICstd.F90
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
       IC file is /glade/u/home/ewinter/derecho/development/kaiju/src/gamera/ICs/null.F90
       Adding executable gamera.x
   Adding Gamera Helio module ...
       IC file is /glade/u/home/ewinter/derecho/development/kaiju/src/gamera/ICs/helio/wsa.F90
       Adding executable gamhelio.x
   Adding ReMIX module ...
       Adding executable remix.x
   Adding RCM module ...
       RCM Grid is of size 180 x 361 x 160
       Adding executable rcm.x
   Adding Voltron module ...
       IC file is /glade/u/home/ewinter/derecho/development/kaiju/src/voltron/ICs/earthcmi.F90
       Adding executable voltron.x
   Adding Base MPI module ...
       Adding executable kaitoy_mpi.x
   Adding Gamera MPI module ...
       Adding executable gamera_mpi.x
   Adding Voltron MPI module ...
       Adding executable voltron_mpi.x
   Adding Gamera Helio MPI module ...
       Adding executable gamhelio_mpi.x
   -- Configuring done (10.9s)
   -- Generating done (1.3s)
   -- Build files have been written to: /glade/u/home/ewinter/derecho/development/kaiju/build_mpi

Step 5: Compile the ``kaiju`` software
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now use ``make`` to build the ``kaiju`` software, time-stamping and saving the output:

.. code-block::

   #!shell
   date; time make -j 4 >& make.out; date


This command should complete in about 10 minutes on ``derecho``. When the command is finished, check the output file ``make.out``. The file is long, but the last few lines should look something like this:

.. code-block::

   [ 98%] Linking Fortran executable bin/voltron_mpi.x
   [ 98%] Built target voltron_mpi.x
   [ 98%] Built target voltron_mpi
   [ 98%] Building Fortran object CMakeFiles/gamhelio_mpi.x.dir/src/gamera/ICs/helio/wsa.F90.o
   [ 99%] Building Fortran object CMakeFiles/gamhelio_mpi.x.dir/src/drivers/gamera_mpix.F90.o
   [100%] Linking Fortran executable bin/gamhelio_mpi.x
   [100%] Built target gamhelio_mpi.x
   [100%] Built target gamhelio_mpi

To verify that all of the ``kaiju`` programs have been built, examine the ``bin`` subdirectory of your ``build_mpi`` directory (this list will evolve as more programs are added):

.. code-block::

   #!shell
   ls bin
   calcdb.x      gamera.x        kaitoy_mpi.x  psd.x   remix2rcm.x    sctrack.x  voltron_mpi.x
   chop.x        gamhelio_mpi.x  kaitoy.x      push.x  remix2remix.x  slice.x    voltron.x
   gamera_mpi.x  gamhelio.x      project.x     rcm.x   remix.x        trace.x    wpicheck.x


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
