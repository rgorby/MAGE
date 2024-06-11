
Building the ``kaiju`` software on a MacOS system
=====================================================

This page provides instructions for building the ``kaiju`` software on a MacOS system. These instructions were developed for MacOS Ventura, but the instructions should work for other MacOS versions.

Introduction
------------

Step 0: Preparing your Mac
--------------------------

A few tools must be installed on your Mac before you begin building the ``kaiju`` software.

Installing HomeBrew
^^^^^^^^^^^^^^^^^^^

`HomeBrew <https://brew.sh/>`_ is a software package management system for MacOS. There are other package management systems available for MacOS. HomeBrew will be used to install tools required for building the ``kaiju`` software.

 If you choose to use an alternate package management system, replace these instructions with the equivalent instructions for your chosen system.

Installing ``cmake``
^^^^^^^^^^^^^^^^^^^^^^^^

`CMake <https://cmake.org/>`_ is used to control the build process for the ``kaiju`` software. This tool is not provided on MacOS systems. You can install it with HomeBrew with the command:

.. code-block::

   #!shell
   brew install cmake


This command will install ``cmake`` and its support files under ``/usr/local`` on your Mac.

Installing Xcode
^^^^^^^^^^^^^^^^

`Xcode <https://developer.apple.com/xcode/>`_ is the standard software development suite for MacOS systems. You can download and install Xcode from the MacOS App Store on your Mac.

Installing the Intel compiler suite
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``kaiju`` software is developed using the Intel compiler suite, available `here <https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html#gs.29x43u>`_. The compiler suite is free, but comes in several optional parts. You will need to install the `Intel oneAPI Base Toolkit <https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html>`_ and the `Intel oneAPI HPC Toolkit <https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html>`_. Once installed, this software should be available on your Mac under ``/opt/intel``.

Once the Intel tools are installed, you will need to make them available during the build process. Add the Intel tools to your command path as follows:

.. code-block::

   #!shell
   export INTEL_HOME=/opt/intel.  # Or wherever you installed it
   export PATH=$INTEL_HOME/oneapi/compiler/latest/mac/bin/intel64:$PATH


A note about directory structures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To ensure compatibility among the ``kaiju`` build instructions for different systems, the instructions for MacOS use a MacOS-specific (currently ``ventura``\ ) subdirectory for all build products. We recommend that you also use this system, to avoid confusion when switching among different computer systems.

Step 1: Building prerequisite libraries
---------------------------------------

NASA CDF library
^^^^^^^^^^^^^^^^

HDF5 library
^^^^^^^^^^^^

Step 2: Creating a python environment
-------------------------------------

Step 3: Compiling the ``kaiju`` software
--------------------------------------------

Serial version
^^^^^^^^^^^^^^

MPI version
^^^^^^^^^^^
