.. role:: raw-html-m2r(raw)
   :format: html


Building the ``kaiju`` software
===================================

These instructions will walk you through the process of building and installing the ``kaiju`` software. Build instructions are provided for MacOS systems, Linux systems, and several HPC (High-Performance Computing) systems (\ ``pleiades``\ , ``cheyenne``\ , and ``frontera``\ ).

**Note:** These instructions presume you have prepared your system and cloned the ``kaiju`` repository, as described `here <./prerequisites.md>`_.

**Note:** Throughout these instructions, ``KAIJUHOME`` refers to the base directory of the `kaiju <https://bitbucket.org/aplkaiju/kaiju>`_ repository that you cloned in the `prerequisites <./prerequisites.md>`_ step.

----

Table of contents
-----------------


* `For all systems <#all>`_
* Personal laptop or desktop machine (serial installation)

  * `MacOS <#macos>`_
  * `Linux <#linux>`_

* HPC serial installation

  * `Cheyenne <#cheyenne-serial>`_
  * `Pleiades <#pleiades-serial>`_

* HPC MPI installation

  * `Cheyenne <#cheyenne-mpi>`_
  * `Pleiades <#pleiades-mpi>`_
  * `Frontera <#frontera-mpi>`_

----

For all systems :raw-html-m2r:`<a name="all"></a>`
------------------------------------------------------

On all systems, you must first create and move to a build directory. The build directory may be anywhere, but it is typically created as a subdirectory named ``build`` under ``KAIJUHOME``\ , like this:

.. code-block::

   #!shell
   mkdir $KAIJUHOME/build
   cd $KAIJUHOME/build


The ``Makefile`` for building ``kaiju`` is created by the ``cmake`` system. When running ``cmake``\ , you may need to include the argument ``-DALLOW_INVALID_COMPILERS=ON`` to force ``cmake`` to use the compiler available on your system.

See `here <tests>`_ for test cases to verify proper installation.
---------------------------------------------------------------------

MacOS (serial) :raw-html-m2r:`<a name="macos"></a>`
-------------------------------------------------------

Intel (x86_64)
--------------

These instructions were developed on MacOS Big Sur, but should work with minimal change on other versions of MacOS.

In your ``build`` directory, run ``cmake`` and build the ``kaiju`` software:

.. code-block::

   #!shell
   FC=`which ifort` cmake -DALLOW_INVALID_COMPILERS=ON $KAIJUHOME
   make
   # Build the GEOPACK library (used by several kaiju scripts).
   cd $KAIJUHOME/external/geopack-2008
   python setup-geopack.py build config_fc --fcompiler=intelem --f77flags=-r8
   # If you use sh/bash/zsh or a compatible shell:
   export PYTHONPATH=$KAIJUHOME/external/geopack-2008/build/lib.macosx-10.9-x86_64-3.8:$PYTHONPATH
   # If you use csh/tcsh or a compatible shell:
   setenv PYTHONPATH $KAIJUHOME/external/geopack-2008/build/lib.macosx-10.9-x86_64-3.8:$PYTHONPATH


**Note:** The path to your ``geopack-2008/build/lib.macosx...`` may differ on your system from the path shown above. Check the contents of the ``build`` directory and use the path appropriate for your system.

**Note:** Monterey: You likely will need to add the Xcode libraries to your LIBRARY_PATH: 

.. code-block::

   #!Shell
   export LIBRARY_PATH="$LIBRARY_PATH:/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib"


If you build HDF5 from source with ifort, you will need to ensure the install path for HDF5 is included in your build settings for kaiju. 

Apple Silicon (M1 Pro, MAX; M2 Pro, MAX)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Using Intel ``ifort`` through Emulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Intel will be deprecating icc in the second-half of 2023 in the oneAPI toolchain. Intel oneAPI DPC++/C++ compiler (\ ``icx``\ ) and oneAPI Fortran Compiler (\ ``ifx``\ ) will not be supported on macOS. You will need to build HDF5 if wanting to use the legacy Intel toolchain (icc and fort) on M1 and newer Macs. Since ifort cannot compile to ``arm64`` you will need to make sure Apple's x86 emulator, `Rosetta 2 <https://support.apple.com/en-us/HT211861>`_\ , is installed. 

.. code-block::

   #!Shell
   softwareupdate --install-rosetta


For an x86_64 build, you will need to be in an emulated terminal mode, the easiest way to ensure this is the case is setting some aliases in you shell config. Example:

.. code-block::

   alias arm="env /usr/bin/arch -arm64 /bin/zsh --login"
   alias intel="env /usr/bin/arch -x86_64 /bin/zsh --login"


To build HDF5 to x86_64 on arm64 Macs it is recommended to use the CMake builds. You will need to change settings in the HDF5options.cmake file with the following (some of these flags are already present, but need to be turned ON or uncommented):

.. code-block::

   #!CMake
   #Turn on Fortran
   set (ADD_BUILD_OPTIONS "${ADD_BUILD_OPTIONS} -DHDF5_BUILD_FORTRAN:BOOL=ON")

   #Select Intel Toolchain
   set (ADD_BUILD_OPTIONS "${ADD_BUILD_OPTIONS} -DCMAKE_TOOLCHAIN_FILE:STRING=config/toolchain/intel.cmake")

   #Set arch to Intel
   set(CMAKE_OSX_ARCHITECTURES "x86_64")


After the build completes, which takes several minutes, you will need to follow the HDF5 CMake install instructions to decide the final installation location and install the HDF5 library and executables. 

If you have h5py installed, make sure you are not in the current venv with that installed. Also, if you have HDF5 installed from Homebrew, ensure that version is not conflicting, or on the PATH, when building ``kaiju`` applications with ``ifort``. 

This HDF5 library and ``h5fc`` wrapper will need to be added in a ``user.cmake`` file in ``$KAIJUHOME/cmake`` to override the default kaiju build. See below for GNU and Intel examples.

.. code-block::

   #!CMake
   #Begin
   if (CMAKE_SYSTEM_NAME MATCHES Darwin)
       # Fix linking on 10.14+. See https://stackoverflow.com/questions/54068035
       LINK_DIRECTORIES(/usr/local/lib)
       SET(CMAKE_OSX_DEPLOYMENT_TARGET 13.0) # Place your Current macOS version number here i.e. 12.0,13.0,14.0 etc. 
   ENDIF()

   if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
       message("User toolchain GNU Selected")
       string(APPEND CMAKE_Fortran_FLAGS " -I /opt/homebrew/Cellar/hdf5/1.14.3/include")
       set(HDF5_LIBRARIES /opt/homebrew/Cellar/hdf5/1.14.3/lib)
       set(HDF5_INCLUDE_DIRS /opt/homebrew/Cellar/hdf5/1.14.3/include)
       set(HDF5_Fortran_COMPILER_EXECUTABLE "/opt/homebrew/Cellar/hdf5/1.14.3/bin/h5fc")
   elseif(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
       set(CMAKE_OSX_ARCHITECTURES $(ARCHS_STANDARD))
       string(APPEND CMAKE_Fortran_FLAGS " -I $USER/HDF5-1.12.2-Darwin/HDF_Group/HDF5/1.12.2/include")
       set(HDF5_LIBRARIES $USER/HDF5-1.12.2-Darwin/HDF_Group/HDF5/1.12.2/lib)
       set(HDF5_INCLUDE_DIRS /$USER/HDF5-1.12.2-Darwin/HDF_Group/HDF5/1.12.2/include)
       set(HDF5_Fortran_COMPILER_EXECUTABLE "$USER/HDF5-1.12.2-Darwin/HDF_Group/HDF5/1.12.2/bin/h5fc")
       set(PKG_CONFIG_PATH $USER/HDF5-1.12.2-Darwin/HDF_Group/HDF5/1.12.2/lib/pkgconfig)
   endif()

   set(CMAKE_Fortran_COMPILER ${HDF5_Fortran_COMPILER_EXECUTABLE})
   set(ALLOW_INVALID_COMPILERS ON)


Depending on your Cmake configuration, the following may be required to be added to your path

.. code-block::

   export PATH="$USER/HDF5-1.12.2-Darwin/HDF_Group/HDF5/1.12.2/bin:$PATH"  



GNU Builds with ``gfortran``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This process is somewhat more straightforward, you will need the above Xcode bits, an install of HDF5 (easiest through Homebrew), and the latest gcc compiler toolchain from Homebrew, or installation from source. If you have added the above cmake options to your ``user.cmake`` file, you should be able to just run the following in your build directory as normal:

.. code-block::

   FC=gfortran cmake ..
   make gamera.x


----

Linux :raw-html-m2r:`<a name="linux"></a>`
----------------------------------------------

TBD

----

Cheyenne (serial) :raw-html-m2r:`<a name="cheyenne-serial"/>`\ </a>
---------------------------------------------------------------------

Environment
^^^^^^^^^^^

Cheyenne uses the ``module`` system to load and unload sets of software from the user environment. Because Cheyenne is sensitive to the order of modules loaded, we suggest purging the current environment modules via 
``module purge`` and then load the latest versions of the required modules. The available versions of each module can be found with the command:

.. code-block::

   #!shell
   module avail


When you load a module, you can specify a specific version, or simply use the latest version of each module. The following modules should be loaded, in the specified order:

.. code-block::

   #!shell
   module purge
   module load git
   module load intel
   module load hdf5   # NOTE: not hdf5-mpi
   module load ncarenv
   module load ncarcompilers
   module load python
   module load cmake


If a specific module version is required, append the version to the module name as ``/version``.

To confirm the currently loaded modules, run ``module list``. For further guidance on modules, see `this reference <https://arc.ucar.edu/knowledge_base/72581272>`_.

Once the modules are loaded, you will need to run Cheyenne's *\ ``ncar_pylib``\ * command to activate the loaded Python module's virtual environment (NPL = NCAR's Python Library). To exit this virtual environment, simply execute *\ ``deactivate``\ *.

If you have created your own environment or want to have multiple collections of modules for various tasks, you can save those `customized environments <https://arc.ucar.edu/knowledge_base/72581272#Environmentmodules-customCustomizedenvironments>`_ for easy re-use.

.. code-block::

   #!shell
   module save environment_name


would save the current module environment for re-use. If we named the environment *kaiju*\ , then we can use ``module restore kaiju`` upon a subsequent login session to restore that working environment.

Compilation
^^^^^^^^^^^

The build system uses ``cmake`` which will attempt to auto-detect HDF5/OpenMP/MPI settings. Optionally, you can provide the file ``cmake/user.cmake`` to set various variables if the auto-detect doesn't work.

.. code-block::

   #!shell
   FC=`which ifort` cmake -DALLOW_INVALID_COMPILERS=ON $KAIJUHOME
   make
   # Build the GEOPACK library (used by several kaiju scripts).
   cd $KAIJUHOME/external/geopack-2008
   python setup-geopack.py build config_fc --fcompiler=intelem --f77flags=-r8
   # If you use sh/bash/zsh or a compatible shell:
   export PYTHONPATH=$KAIJUHOME/external/geopack-2008/build/lib.linux-x86_64-3.8:$PYTHONPATH
   # If you use csh/tcsh or a compatible shell:
   setenv PYTHONPATH $KAIJUHOME/external/geopack-2008/build/lib.linux-x86_64-3.8:$PYTHONPATH


**Note:** The path to your ``geopack-2008/build/lib...`` may differ on your system from the path shown above. Check the contents of the ``build`` directory and use the path appropriate for your system.

 **Note:** Running ``make`` or ``make all`` will build all available components and configurations and can take a considerable amount of time to build and file storage. You can also specify individual components to build. To get the complete list of build targets you can issue the ``make help`` command.

The ``bin`` subdirectory of the build directory now contains many different executables, such as ``gamera.x``. See `here <tests>`_ for test cases to verify proper installation.

----

Pleiades (serial) :raw-html-m2r:`<a name="pleiades-serial"></a>`
--------------------------------------------------------------------

Loading required environment modules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Pleiades uses the ``module`` system to load and unload sets of software from the user environment. Because Pleiades is sensitive to the order of modules loaded, we suggest purging the current environment modules via ``module purge`` and then load the latest versions of the required modules. The available versions of each module can be found with the command:

.. code-block::

   #!shell
   module avail


When you load a module, you can specify a specific version, or simply use the latest version of each module. For Pleiades, we suggest the following modules should be loaded, in the specified order:

.. code-block::

   #!shell
   module purge
   module load pkgsrc/2021Q2
   module load comp-intel/2020.4.304
   module load szip/2.1.1
   module load hdf5/1.8.18_serial


These exact modules must be used (as of 3 March 2022).

Building the code
^^^^^^^^^^^^^^^^^

The ``cmake`` process may not correctly identify the Intel Fortran compiler (\ ``ifort``\ ) as the desired Fortran compiler, so the ``FC`` environment variable must be specified on the ``cmake`` command line in order to use the correct compiler. The following steps perform these actions.

.. code-block::

   #!shell
   # The cmake command takes a few seconds.
   FC=`which ifort` cmake $KAIJUHOME
   # The make command can take 20 minutes or more on pleiades.
   make
   # Build the GEOPACK library (used by several kaiju scripts).
   cd $KAIJUHOME/external/geopack-2008
   python setup-geopack.py build config_fc --fcompiler=intelem --f77flags=-r8
   # If you use sh/bash/zsh or a compatible shell:
   export PYTHONPATH=$KAIJUHOME/external/geopack-2008/build/lib.linux-x86_64-3.8:$PYTHONPATH
   # If you use csh/tcsh or a compatible shell:
   setenv PYTHONPATH $KAIJUHOME/external/geopack-2008/build/lib.linux-x86_64-3.8:$PYTHONPATH


**Note:** The path to your ``geopack-2008/build/lib...`` may differ on your system from the path shown above. Check the contents of the ``build`` directory and use the path appropriate for your system.

When the build is complete, you will find the compiled executables in the ``bin`` subdirectory of your build directory:

.. code-block::

   #!shell
   bash-4.2$ ls bin
   calcdb.x  gamhelio.x  psd.x   remix.x        sctrack.x  voltron.x
   chop.x    kaitoy.x    push.x  remix2rcm.x    slice.x    wpicheck.x
   gamera.x  project.x   rcm.x   remix2remix.x  trace.x


You can add this directory to your ``PATH``\ :

.. code-block::

   #!shell
   # If you use sh/bash/zsh or a compatible shell:
   export PATH=$KAIJUHOME/build/bin:$PATH
   # If you use csh/tcsh or a compatible shell:
   setenv PATH $KAIJUHOME/build/bin:$PATH



See `here <tests>`_ for test cases to verify proper installation.

----

Cheyenne (MPI) :raw-html-m2r:`<a name="cheyenne-mpi"></a>`
--------------------------------------------------------------

Environment
^^^^^^^^^^^

Cheyenne uses the ``module`` system to load and unload sets of software from the user environment. Because Cheyenne is sensitive to the order of modules loaded, we suggest purging the current environment modules via ``module purge`` and then load the latest versions of the required modules. The available versions of each module can be found with the command:

.. code-block::

   #!shell
   module avail


When you load a module, you can specify a specific version, or simply use the latest version of each module. The following modules should be loaded, in the specified order:

.. code-block::

   #!shell
   module purge
   module load git
   module load intel
   module load impi
   module load hdf5-mpi
   module load ncarenv
   module load ncarcompilers
   module load python
   module load cmake


To confirm the currently loaded modules, run ``module list``. For further guidance on modules, see `this reference <https://arc.ucar.edu/knowledge_base/72581272>`_.

Once the modules are loaded, you will need to run Cheyenne's *\ ``ncar_pylib``\ * command to activate the loaded Python module's virtual environment (NPL = NCAR's Python Library). To exit this virtual environment, simply execute *\ ``deactivate``\ *.

Building the code
^^^^^^^^^^^^^^^^^

In your build directory, run the following commands:

.. code-block::

   #!shell
   # The cmake command takes a few seconds.
   cmake -DENABLE_MPI=ON $KAIJUHOME
   # The make command can take 20 minutes or more on cheyenne.
   make
   # Build the GEOPACK library (used by several kaiju scripts).
   cd $KAIJUHOME/external/geopack-2008
   python setup-geopack.py build config_fc --fcompiler=intelem --f77flags=-r8
   # If you use sh/bash/zsh or a compatible shell:
   export PYTHONPATH=$KAIJUHOME/external/geopack-2008/build/lib.linux-x86_64-3.8:$PYTHONPATH
   # If you use csh/tcsh or a compatible shell:
   setenv PYTHONPATH $KAIJUHOME/external/geopack-2008/build/lib.linux-x86_64-3.8:$PYTHONPATH


When the build is complete, you will find the compiled executables in the bin subdirectory of your build directory:

.. code-block::

   #!shell
   bash-4.2$ ls bin
   calcdb.x  kaitoy.x   push.x   remix2rcm.x    slice.x    wpicheck.x
   chop.x    project.x  rcm.x    remix2remix.x  trace.x
   gamera.x  psd.x      remix.x  sctrack.x      voltron.x


Add this directory to your ``PATH``\ :

.. code-block::

   #!shell
   # If you use sh/bash/zsh or a compatible shell:
   export PATH=$KAIJUHOME/build/bin:$PATH
   # If you use csh/tcsh or a compatible shell:
   setenv PATH $KAIJUHOME/build/bin:$PATH


See `here <tests>`_ for test cases to verify proper installation.

Pleiades (MPI) :raw-html-m2r:`<a name="pleiades-mpi"></a>`
--------------------------------------------------------------

Loading required environment modules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Pleiades uses the ``module`` system to load and unload sets of software from the user environment. Because Pleiades is sensitive to the order of modules loaded, we suggest purging the current environment modules via ``module purge`` and then load the latest versions of the required modules. The available versions of each module can be found with the command:

.. code-block::

   #!shell
   module avail


When you load a module, you can specify a specific version, or simply use the latest version of each module. For Pleiades, we suggest the following modules should be loaded, in the specified order:

.. code-block::

   #!shell
   module purge
   module load pkgsrc/2021Q2
   module load comp-intel/2020.4.304
   module load mpi-hpe/mpt.2.23
   module load hdf5/1.8.18_mpt


These exact modules must be used (as of 3 March 2022), since the MPI libraries are built specifically for the Intel compilers - MPI for ``gcc`` is not available yet on pleiades.

Building the code
^^^^^^^^^^^^^^^^^

The ``cmake`` process may not correctly identify the Intel fortran compiler (\ ``ifort``\ ) as the desired Fortran compiler, so the ``FC`` environment variable must be specified on the ``cmake`` command line in order to use the correct compiler. The following steps perform these actions.

.. code-block::

   #!shell
   # The cmake command takes a few seconds.
   FC=ifort cmake -DENABLE_MPI=ON ..
   # The make command can take 20 minutes or more on pleiades.
   make
   # Build the GEOPACK library (used by several kaiju scripts).
   cd $KAIJUHOME/external/geopack-2008
   python setup-geopack.py build config_fc --fcompiler=intelem --f77flags=-r8
   # If you use sh/bash/zsh or a compatible shell:
   export PYTHONPATH=$KAIJUHOME/external/geopack-2008/build/lib.linux-x86_64-3.8:$PYTHONPATH
   # If you use csh/tcsh or a compatible shell:
   setenv PYTHONPATH $KAIJUHOME/external/geopack-2008/build/lib.linux-x86_64-3.8:$PYTHONPATH


When the build is complete, you will find the compiled executables in the ``bin`` subdirectory of your build directory:

.. code-block::

   #!shell
   bash-4.2$ ls bin
   calcdb.x      gamhelio_mpi.x  push.x         sctrack.x      wpicheck.x
   chop.x        kaitoy.x        rcm.x          slice.x
   gamera.x      kaitoy_mpi.x    remix.x        trace.x
   gamera_mpi.x  project.x       remix2rcm.x    voltron.x
   gamhelio.x    psd.x           remix2remix.x  voltron_mpi.x


Add this directory to your ``PATH``\ :

.. code-block::

   #!shell
   # If you use sh/bash/zsh or a compatible shell:
   export PATH=$KAIJUHOME/build/bin:$PATH
   # If you use csh/tcsh/zsh or a compatible shell:
   setenv PATH $KAIJUHOME/build/bin:$PATH


See `here <tests>`_ for test cases to verify proper installation.

Frontera (MPI) :raw-html-m2r:`<a name="frontera-mpi"></a>`
--------------------------------------------------------------

Loading required environment modules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Execute the following ``module`` commands to load the various compilers and libraries needed to build the ``kaiju`` code on Frontera:

.. code-block::

   #!shell
   module load intel/19.1.3
   module load phdf5


This will give a module environment that includes:

.. code-block::

   #!shell
   1) git/2.24.1      4) pmix/3.1.4      7) TACC          10) python3/3.7.0
   2) autotools/1.2   5) hwloc/1.11.12   8) intel/19.1.3  11) phdf5/1.10.4
   3) cmake/3.20.3    6) xalt/2.10.31    9) impi/19.0.9


**???** As of the writing (11/2021) the ``intel/19.1.3`` module was not publicly available on Frontera so you will have to download it. 

Building the code
^^^^^^^^^^^^^^^^^

In the folder ``$KAIJUHOME/cmake`` create a file "user.cmake" to set some compiler and linker flags as below:

.. code-block::

   #!shell
   string(APPEND CMAKE_Fortran_FLAGS " -I$ENV{TACC_HDF5_INC}")
   string(APPEND CMAKE_Fortran_FLAGS " -L$ENV{TACC_HDF5_LIB} -lhdf5_fortran -lhdf5hl_fortran")


Once this step is complete, you can build the ``kaiju`` software using the Intel compilers. The ``cmake`` process does not correctly identify the Intel fortran compiler (\ ``ifort``\ ) as the desired Fortran compiler, so the ``FC`` environment variable must be specified on the ``cmake`` command line in order to use the correct compiler.

.. code-block::

   #!shell
   # The cmake command takes a few seconds.
   FC=h5pfc cmake -DENABLE_MPI=ON $KAIJUHOME
   make
   # Build the GEOPACK library (used by several kaiju scripts).
   cd $KAIJUHOME/external/geopack-2008
   python setup-geopack.py build config_fc --fcompiler=intelem --f77flags=-r8
   # If you use sh/bash/zsh or a compatible shell:
   export PYTHONPATH=$KAIJUHOME/external/geopack-2008/build/lib.linux-x86_64-3.8:$PYTHONPATH
   # If you use csh/tcsh or a compatible shell:
   setenv PYTHONPATH $KAIJUHOME/external/geopack-2008/build/lib.linux-x86_64-3.8:$PYTHONPATH


When the build is complete, you will find the compiled executables in the bin subdirectory of your build directory:

.. code-block::

   #!shell
   bash-4.2$ ls bin
   calcdb.x  kaitoy.x   push.x   remix2rcm.x    slice.x    wpicheck.x
   chop.x    project.x  rcm.x    remix2remix.x  trace.x
   gamera.x  psd.x      remix.x  sctrack.x      voltron.x


Add this directory to your ``PATH``\ :

.. code-block::

   #!shell
   # If you use sh/bash/zsh or a compatible shell:
   export PATH=$KAIJUHOME/build/bin:$PATH
   # If you use csh/tcsh or a compatible shell:
   setenv PATH $KAIJUHOME/build/bin:$PATH


**Note:** When using the MPT, there are additional scripts and environment variables that must be set. Please consult the `Running with MPT <runningWithMPT>`_ wiki page for information on using MPT. If using Intel MPI, these additional steps are not necessary.

See `here <tests>`_ for test cases to verify proper installation.
