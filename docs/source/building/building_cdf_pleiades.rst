
Introduction
------------

This file describes how to build and install the CDF library on ``pleiades`` for use with the ``kaiju`` software.

As of this date (2023-06-30 12:24:55), the latest version of the CDF library is:

``cdf39_0-dist-all.tar.gz``

NOTE: A ``pleiades``\ -specific subdirectory is used.

.. code-block::

   #!shell
   # NOTE: No modules are loaded on pleiades by default.

   # Specify the name for this machine.
   export HOST_SYSTEM=pleiades

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
   # Took 0m46.197s

   # Test the library.
   date; time make test >& make_test.out; date
   # Took 0m0.240s

   # Install the library in a version-specific subdirectory.
   date; time make INSTALLDIR=$BUILD_ROOT install >& make_install.out; date
   # Took 0m1.413s

   # Clean the build tree.
   date; time make clean >& make_clean.out; date
   # Took 0m0.226s


To use this software, you must run the setup script:

.. code-block::

   #!shell
   source $BUILD_ROOT/bin/definitions.B
