
Building and installing the NASA CDF library on ``pleiades``
================================================================

Introduction
------------

This file describes how to build and install the CDF library on ``pleiades`` for use with the ``kaiju`` software.

These instructions describe use of CDF 3.9.0. The software will be built and installed into:

.. code-block:: shell

   $HOME/local/cdf/3.9.0


Building and installing the CDF library
---------------------------------------

.. code-block:: shell

   # NOTE: No modules are loaded on pleiades by default.

   # Specify and create the root of the build tree.
   export BUILD_ROOT=$HOME/local/cdf/3.9.0
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
   # Took 1m54.712s

   # Test the library.
   date; time make test >& make_test.out; date
   # Took 0m3.563s

   # Install the library in a version-specific subdirectory.
   date; time make INSTALLDIR=$BUILD_ROOT install >& make_install.out; date
   # Took 0m13.776s

   # Clean the build tree.
   date; time make clean >& make_clean.out; date
   # Took 0m2.718s


Using the CDF library
=====================

To use this software, you must run the setup script:

.. code-block:: shell

   source $BUILD_ROOT/bin/definitions.B
