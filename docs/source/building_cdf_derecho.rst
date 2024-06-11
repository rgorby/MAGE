
Introduction
------------

This file describes how to build and install the CDF library on ``derecho`` for use with the ``kaiju`` software.

*NOTE*\ : A ``cheyenne``\ -specific subdirectory of the user home directory is used for this build.

As of this date (2023-06-27 15:48:52), the latest version of the CDF library is:

``cdf39_0-dist-all.tar.gz``

.. code-block::

   #!shell
   # USING THE DEFAULT MODULE SET ON derecho, WHICH IS:
   # 1) ncarenv/23.06
   # 2) craype/2.7.20
   # 3) intel/2023.0.0
   # 4) ncarcompilers/1.0.0
   # 5) cray-mpich/8.1.25
   # 6) hdf5/1.12.2
   # 7) netcdf/4.9.2

   # Specify the name for this machine.
   export HOST_SYSTEM=derecho

   # Specify and create the root of the build tree.
   export BUILD_ROOT=$HOME/$HOST_SYSTEM/local/cdf/3.9.0
   mkdir -p $BUILD_ROOT/src
   cd $BUILD_ROOT/src

   # Download the source tarball.
   wget https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/cdf39_0/linux/cdf39_0-dist-all.tar.gz

   # Unpack the source code.
   tar xzvf cdf39_0-dist-all.tar.gz

   # Move into the code directory.
   cd cdf39_0-dist

   # Build the library using the default system GNU compiler.
   # NOTE: curses not available in default derecho environment.
   date; time make OS=linux ENV=gnu CURSES=no all >& make.out; date
   # Took 0m45.548s

   # Test the library.
   date; time make test >& make_test.out; date
   # Took 0m0.423s

   # Install the library in a version-specific subdirectory.
   date; time make INSTALLDIR=$BUILD_ROOT install >& make_install.out; date
   # Took 0m2.261s

   # Clean the build tree.
   date; time make clean >& make_clean.out; date
   # Took 0m0.394s


To use this software, you must run the setup script:

.. code-block::

   #!shell
   source $BUILD_ROOT/bin/definitions.B
