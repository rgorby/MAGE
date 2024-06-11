
This file describes how to build and install the CDF library on ``cheyenne`` for use with the ``kaiju`` software.

*NOTE*\ : A ``cheyenne``\ -specific subdirectory of the user home directory is used for this build.

As of this date (2023-06-27 15:48:52), the latest version of the CDF library is:

``cdf39_0-dist-all.tar.gz``

.. code-block::

   #!shell
   # USING THE DEFAULT MODULE SET ON CHEYENNE, WHICH IS:
   # 1) ncarenv/1.3
   # 2) intel/19.1.1
   # 3) ncarcompilers/0.5.0
   # 4) mpt/2.25
   # 5) netcdf/4.8.1
   # 6) tmux/2.9a   # Added by me

   # Specify the name for this machine.
   export HOST_SYSTEM=cheyenne

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
   date; time make OS=linux ENV=gnu all >& make.out; date
   # Took 0m50.743s

   # Test the library.
   date; time make test >& make_test.out; date
   # Took 0m0.726s

   # Install the library in a version-specific subdirectory.
   date; time make INSTALLDIR=$BUILD_ROOT install >& make_install.out; date
   # Took 0m3.553s

   # Clean the build tree.
   date; time make clean >& make_clean.out; date
   # Took 0m0.414s


To use this software, you must run the setup script:

.. code-block::

   #!shell
   source $BUILD_ROOT/bin/definitions.B
