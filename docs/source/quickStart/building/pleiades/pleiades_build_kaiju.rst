
Building the ``kaiju`` software on the ``pleiades`` HPC system
======================================================================

Introduction
------------

These instructions will walk you through the process of building and installing the ``kaiju`` software on the ``pleiades`` supercomputer.

These instructions assume that the user is using the ``bash`` shell, and that no modifications have been made to the user "dotfiles" (\ ``$HOME/.bashrc``\ , ``$HOME/.bash_profile``\ ). If you have customized either of these files for your account, please carefully save and inspect the output from each command in the build process to ensure that no unexpected problems have crept in. To facilitate this practice, all of the commands shown below will illustrate how to save command output, and how to measure how long each step takes. The latter is a useful bit of information which can help identify build problems early in the process, avoiding much wasted time and effort later.

Like most HPC systems, ``pleiades`` uses the ``module`` system to manage the versions of software packages available to the user. A "module" is a collection of programs and libraries for a specific task, and "loading" the module adjusts the user environment (mostly by setting or updating environment variables) to make that module available to the user. When you log in to ``pleiades``\ , no modules are loaded by default (NOTE: This is new behavior in pleiades since the TOSS4 upgrade.):

.. code-block:: shell

   module list
   No Modulefiles Currently Loaded.


The correct module sets for the serial and MPI versions of ``kaiju`` will be listed in the instructions below.

Building the serial version of the ``kaiju`` software on ``pleiades``
-----------------------------------------------------------------------------

Step 1: Load the build modules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Begin by loading the modules needed to build the serial version of the ``kaiju`` software on ``pleiades``\ :

.. code-block:: shell

   module purge
   module load nas
   module load pkgsrc/2022Q1-rome  # For git-lfs and cmake
   module load comp-intel/2020.4.304  # Latest version
   module load szip/2.1.1
   module load hdf5/1.8.18_serial


After these commands have been run, verify your module list with ``module list``\ :

.. code-block:: shell

   module list

   Currently Loaded Modules:
    1) nas   2) pkgsrc/2022Q1-rome   3) comp-intel/2020.4.304   4) szip/2.1.1   5) hdf5/1.8.18_serial  


NOTE: This list of modules will evolve. The current list is the recommended module set as of 2023-06-30.

Step 2: Create the build directory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this case, just move to your home directory.

.. code-block:: shell

   cd $HOME


Step 3: Clone the ``kaiju`` repository from BitBucket
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

NOTE: This step assumes you have been granted access to the ``kaiju`` repository on BitBucket, and that you have configured an SSH key pair for use with BitBucket. If you need help with these tasks, please contact a CGS team member for assistance.

Clone the ``kaiju`` repository (or "repo") from BitBucket:

.. code-block:: shell

   git clone git@bitbucket.org:aplkaiju/kaiju.git


This process should take a minute or so. When complete, verify that the ``kaiju`` code exists in your directory (the actual directory contents may differ slightly from what is shown below):

.. code-block:: shell

   ls kaiju
   analysis        examples        kaiju.sublime-project  pytests     scripts   testingScripts
   cmake           external        kaipy                  quickstart  setup.py  tests
   CMakeLists.txt  gitHookScripts  places                 README.md   src       xml


Now move down into the cloned repo, and switch to the branch of the code you wish to use. By default, the cloned repository provides the ``master`` branch, but we want the ``development`` branch:

.. code-block:: shell

   cd kaiju
   git switch development


Step 4: Run ``cmake`` to create the ``Makefile`` needed to build the software
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Since the ``kaiju`` code can be built in serial and MPI forms, we first make a directory in which to build the serial version of the code (use whatever name you prefer, but ``build_serial`` is simple and unambiguous):

.. code-block:: shell

   mkdir build_serial
   cd build_serial


Now run the ``cmake`` command. Save the ``cmake`` output, and use timestamps for each step. The options shown below direct the build process to use a recent version of the Intel Fortran compiler:

.. code-block:: shell

   date; time FC=`which ifort` cmake >& cmake.out; date


This command usually takes about 5 seconds, depending on system activity. Examine the output file ``cmake.out`` for problems. It *should* look something like this:

.. code-block::

   -- The Fortran compiler identification is Intel 19.1.3.20200925
   -- Detecting Fortran compiler ABI info
   -- Detecting Fortran compiler ABI info - done
   -- Check for working Fortran compiler: /nasa/intel/Compiler/2020.4.304/compilers_and_libraries_2020.4.304/linux/bin/intel64/ifort - skipped
   -- Found HDF5: /nasa/hdf5/1.8.18_serial/lib/libhdf5_fortran.so;/nasa/hdf5/1.8.18_serial/lib/libhdf5.so;/nasa/szip/2.1.1/lib/libsz.so;/nasa/pkgsrc/toss4/2022Q1-rome/lib/libz.so;/usr/lib64/libdl.so;/usr/lib64/libm.so (found version "1.8.18") found componen\
   ts: Fortran
   -- Found OpenMP_Fortran: -qopenmp (found version "5.0")
   -- Found OpenMP: TRUE (found version "5.0") found components: Fortran
   -------------------------
   Configuration summary ...
   System: pfe26
   OS: Linux
   Processor: x86_64
   Compiler: Intel / 19.1.3.20200925
             /nasa/intel/Compiler/2020.4.304/compilers_and_libraries_2020.4.304/linux/bin/intel64/ifort
   HDF5 Wrapper: /nasa/hdf5/1.8.18_serial/bin/h5fc
   Version: 3f4e147c / development
   Build Type: Release
   Base Flags:  -fPIC -free -implicitnone -qopenmp
   Build Flags: -O3 -align array64byte -align rec32byte -no-prec-div -fast-transcendentals -ipo -march=corei7 -axCORE-AVX2
   -------------------------

   Adding CHIMP module ...
       EB IC file is /home3/ewinter/pleiades/cgs/kaiju/development/kaiju/src/chimp/ebICs/ebICstd.F90
       TP IC file is /home3/ewinter/pleiades/cgs/kaiju/development/kaiju/src/chimp/tpICs/tpICstd.F90
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
       IC file is /home3/ewinter/pleiades/cgs/kaiju/development/kaiju/src/gamera/ICs/null.F90
       Adding executable gamera.x
   Adding Gamera Helio module ...
       IC file is /home3/ewinter/pleiades/cgs/kaiju/development/kaiju/src/gamera/ICs/helio/wsa.F90
       Adding executable gamhelio.x
   Adding ReMIX module ...
       Adding executable remix.x
   Adding RCM module ...
       RCM Grid is of size 180 x 361 x 160
       Adding executable rcm.x
   Adding Voltron module ...
       IC file is /home3/ewinter/pleiades/cgs/kaiju/development/kaiju/src/voltron/ICs/earthcmi.F90
       Adding executable voltron.x
   -- Configuring done
   -- Generating done
   -- Build files have been written to: /home3/ewinter/pleiades/cgs/kaiju/development/kaiju/build_serial

Step 5: Compile the ``kaiju`` software
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now use ``make`` to build the ``kaiju`` software, time-stamping and saving the output:

.. code-block:: shell

   date; time make >& make.out; date


This command should complete in about 20 minutes. When the command is finished, check the output file ``make.out``. The file is long, but the last few lines should look something like this:

.. code-block::

   /home3/ewinter/pleiades/cgs/kaiju/development/kaiju/src/voltron/ICs/earthcmi.F90(48): remark #15009: uservoltic_mp_inituser_ has been targeted for automatic cpu dispatch
   /home3/ewinter/pleiades/cgs/kaiju/development/kaiju/src/voltron/ICs/earthcmi.F90(169): remark #15009: uservolticinituser_mp_psphereic_ has been targeted for automatic cpu dispatch
   /home3/ewinter/pleiades/cgs/kaiju/development/kaiju/src/voltron/voltapp.F90(28): remark #15009: voltapp_mp_initvoltron_ has been targeted for automatic cpu dispatch
   /home3/ewinter/pleiades/cgs/kaiju/development/kaiju/src/base/kronos.F90(39): remark #15009: kronos_mp_initts_ has been targeted for automatic cpu dispatch
   [100%] Built target voltron.x
   [100%] Built target voltron

To verify that all of the ``kaiju`` programs have been built, examine the ``bin`` subdirectory of your ``build_serial`` directory (this list will evolve as more programs are added):

.. code-block:: shell

   ls bin
   calcdb.x  chop.x  gamera.x  gamhelio.x  kaitoy.x  project.x  psd.x  push.x  rcm.x  remix2rcm.x  remix2remix.x  remix.x  sctrack.x  slice.x  trace.x  voltron.x  wpicheck.x




Building the MPI version of the ``kaiju`` software on ``pleiades``
--------------------------------------------------------------------------

Step 1: Load the build modules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Begin by loading the modules needed to build the MPI version of the ``kaiju`` software on ``pleiades``\ :

.. code-block:: shell

   module load nas
   module load pkgsrc/2022Q1-rome  # For git-lfs and cmake
   module load comp-intel/2020.4.304  # Latest version
   module load mpi-hpe/mpt.2.23
   module load szip/2.1.1
   module load hdf5/1.8.18_mpt


After these commands have been run, verify your module list with ``module list``\ :

.. code-block:: shell

   module list

   Currently Loaded Modules:
    1) nas   2) pkgsrc/2022Q1-rome   3) comp-intel/2020.4.304   4) mpi-hpe/mpt.2.23   5) szip/2.1.1   6) hdf5/1.8.18_mpt  


NOTE: This list of modules will evolve. The current list is the recommended module set as of 2023-06-30.

Step 2: Create the build directory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this case, just move to your home directory.

.. code-block:: shell

   cd $HOME


Step 3: Clone the ``kaiju`` repository from BitBucket
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

NOTE: This step assumes you have been granted access to the ``kaiju`` repository on BitBucket, and that you have configured an SSH key pair for use with BitBucket. If you need help with these tasks, please contact a CGS team member for assistance.

Clone the ``kaiju`` repository (or "repo") from BitBucket (skip this part if you already cloned the repo):

.. code-block:: shell

   git clone git@bitbucket.org:aplkaiju/kaiju.git


This process should take a minute or so. When complete, verify that the ``kaiju`` code exists in your directory (the actual directory contents may differ slightly from what is shown below):

.. code-block:: shell

   ls kaiju
   analysis      CMakeLists.txt  gitHookScripts         places      README.md  src             xml
   build_serial  examples        kaiju.sublime-project  pytests     scripts    testingScripts
   cmake         external        kaipy                  quickstart  setup.py   tests


Now move down into the cloned repo, and switch to the branch of the code you wish to use. By default, the cloned repository provides the ``master`` branch, but we want the ``development`` branch:

.. code-block:: shell

   cd kaiju
   git switch development


Step 4: Run ``cmake`` to create the ``Makefile`` needed to build the software
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Since the ``kaiju`` code can be built in serial and MPI forms, we first make a directory in which to build the MPI version of the code (use whatever name you prefer, but ``build_mpi`` is simple and unambiguous):

.. code-block:: shell

   mkdir build_mpi
   cd build_mpi


Now run the ``cmake`` command. Save the ``cmake`` output, and use timestamps for each step. The options shown below direct the build process to use a recent version of the Intel Fortran compiler:

.. code-block:: shell

   date; time FC=`which ifort` cmake -DENABLE_MPI=ON .. >& cmake.out; date


This command usually takes 5-10 seconds, depending on system activity. Examine the output file ``cmake.out`` for problems. It *should* look something like this:

.. code-block::

   Using default MPI Address Size of MPI_ADDRESS_KIND for basic MPI functions
   Using default MPI Address Size of MPI_ADDRESS_KIND for neighborhood MPI functions
   CMake Warning at /nasa/pkgsrc/toss4/2022Q1-rome/share/cmake-3.22/Modules/FindHDF5.cmake:745 (message):
     HDF5 found for language Fortran is parallel but previously found language
     is not parallel.
   Call Stack (most recent call first):
     cmake/compilers.cmake:5 (find_package)
     CMakeLists.txt:99 (include)


   -- Found HDF5: /nasa/hdf5/1.8.18_mpt/lib/libhdf5_fortran.so;/nasa/hdf5/1.8.18_mpt/lib/libhdf5.so;/nasa/szip/2.1.1/lib/libsz.so;/nasa/pkgsrc/toss4/2022Q1-rome/lib/libz.so;/usr/lib64/libdl.so;/usr/lib64/libm.so (found version "1.8.18") found components: Fo\
   rtran
   -- Found OpenMP_Fortran: -qopenmp (found version "5.0")
   -- Found OpenMP: TRUE (found version "5.0") found components: Fortran
   -- Found MPI
   -------------------------
   Configuration summary ...
   System: pfe26
   OS: Linux
   Processor: x86_64
   Compiler: Intel / 19.1.3.20200925
             mpif08
   HDF5 Wrapper: /nasa/hdf5/1.8.18_mpt/bin/h5pfc
   Version: 3f4e147c / development
   Build Type: Release
   Base Flags:  -fPIC -free -implicitnone -qopenmp  -lmpi
   Build Flags: -O3 -align array64byte -align rec32byte -no-prec-div -fast-transcendentals -ipo -march=corei7 -axCORE-AVX2
   -------------------------

   Adding CHIMP module ...
       EB IC file is /home3/ewinter/pleiades/cgs/kaiju/development/kaiju/src/chimp/ebICs/ebICstd.F90
       TP IC file is /home3/ewinter/pleiades/cgs/kaiju/development/kaiju/src/chimp/tpICs/tpICstd.F90
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
       IC file is /home3/ewinter/pleiades/cgs/kaiju/development/kaiju/src/gamera/ICs/null.F90
       Adding executable gamera.x
   Adding Gamera Helio module ...
       IC file is /home3/ewinter/pleiades/cgs/kaiju/development/kaiju/src/gamera/ICs/helio/wsa.F90
       Adding executable gamhelio.x
   Adding ReMIX module ...
       Adding executable remix.x
   Adding RCM module ...
       RCM Grid is of size 180 x 361 x 160
       Adding executable rcm.x
   Adding Voltron module ...
       IC file is /home3/ewinter/pleiades/cgs/kaiju/development/kaiju/src/voltron/ICs/earthcmi.F90
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
   -- Build files have been written to: /home3/ewinter/pleiades/cgs/kaiju/development/kaiju/build_mpi

Step 5: Compile the ``kaiju`` software
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now use ``make`` to build the ``kaiju`` software, time-stamping and saving the output:

.. code-block:: shell

   date; time make >& make.out; date


This command should complete in about 30 minutes. When the command is finished, check the output file ``make.out``. The file is long, but the last few lines should look something like this:

.. code-block::

   /home3/ewinter/pleiades/cgs/kaiju/development/kaiju/src/gamera/ICs/helio/wsa.F90(70): remark #15009: usergamic_mp_inituser_ has been targeted for automatic cpu dispatch                                                                                      
   /home3/ewinter/pleiades/cgs/kaiju/development/kaiju/src/gamera/ICs/helio/wsa.F90(222): remark #15009: usergamicinituser_mp_mapkinit_ has been targeted for automatic cpu dispatch                                                                             
   /home3/ewinter/pleiades/cgs/kaiju/development/kaiju/src/gamera/ICs/helio/wsa.F90(201): remark #15009: usergamicinituser_mp_gasic_ has been targeted for automatic cpu dispatch                                                                                
   [100%] Built target gamhelio_mpi.x
   [100%] Built target gamhelio_mpi

To verify that all of the ``kaiju`` programs have been built, examine the ``bin`` subdirectory of your ``build_mpi`` directory (this list will evolve as more programs are added):

.. code-block:: shell

   ls bin
   calcdb.x  chop.x  gamera_mpi.x  gamera.x  gamhelio_mpi.x  gamhelio.x  kaitoy_mpi.x  kaitoy.x  project.x  psd.x  push.x  rcm.x  remix2rcm.x  remix2remix.x  remix.x  sctrack.x  slice.x  trace.x  voltron_mpi.x  voltron.x  wpicheck.x


Using the ``kaiju`` software
--------------------------------

Once built, you must run the setup script before using the ``kaiju`` software:

.. code-block:: shell

   source $KAIJU_HOME/scripts/setupEnvironment.sh


This script will set environment variables needed by the ``kaiju`` software, including the ``KAIJUHOME`` environment variable (not the ``KAIJU_HOME`` environment variable). However, the path to the compiled programs is not added - you will need to specify the complete path when using compiled programs. For example,  to run the serial version of ``gamera.x``\ :

.. code-block:: shell

   $KAIJUHOME/build_serial/bin/gamera.x
