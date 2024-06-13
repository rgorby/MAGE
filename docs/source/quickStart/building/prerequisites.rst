
Prerequisites for building the kaiju software
=============================================

----

Introduction
------------

The ``kaiju`` software is a combination of Fortran and Python code. In general, Python is used for pre- and post-processing, data analysis, and various utility functions, while Fortran executables are used as the computational core. The ``kaiju`` software was designed to run on individual laptop computers as well as HPC systems.

The ``kaiju`` code is managed on BitBucket at https://bitbucket.org/aplkaiju/kaiju.

This page describes the prerequisites that your computer system must provide in order to build and run the ``kaiju`` software. These instructions assume you are comfortable using a command-line interface: a command shell (\ ``bash`` is assumed), compilers, Python, git, and other standard tools.

:pushpin: **Note**\ : Throughout these instructions, ``KAIJU_HOME`` refers to the base directory of the ``kaiju`` code, created when you cloned the `\ ``kaiju`` repository <https://bitbucket.org/aplkaiju/kaiju/src/master/>`_ (see below). For example, if you cloned the ``kaiju`` repository in your home directory, e.g. ``/home/username``\ , then ``KAIJU_HOME`` would be ``/home/username/kaiju``.

----

Fortran compiler
----------------

The ``kaiju`` code is developed primarily using the Intel Fortran compiler, version 17 or later, which is available as a `free download from Intel <https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html#gs.3msomg>`_. The code can also be compiled using the GNU Fortran compiler.

If you need to install the Intel Fortran compiler, instructions can be found :doc:`here <install_intel_compilers>`.

----

NASA CDF library
----------------

The SpacePy package used by the Python portion of the ``kaiju`` software requires the `NASA CDF library <https://cdf.gsfc.nasa.gov/>`_. If you need to install the CDF library, instructions can be found :doc:`here </build_cdf>`.

----

GEOS library
------------

The Python portion of the ``kaiju`` software uses the CartoPy package, which requires the `GEOS library <https://libgeos.org/>`_. This library is typically available as an environment module on HPC systems. However, the GEOS library is not available on ``pleiades``\ , so ``kaiju`` code which uses GEOS will not work on ``pleiades``. On a Mac, the GEOS library can be installed using HomeBrew:

.. code-block::

   #!shell
   brew install geos



----

HDF5 library
------------

The ``kaiju`` software makes extensive use of the `HDF5 file format <https://www.hdfgroup.org/>`_ for data files, and so requires an HDF5 library for proper operation. Many systems provide the HDF5 library by default, or via a system-wide mechanism (such as the ``module`` system on many HPC systems).

If the HDF5 library is not available, you must install it. The code can usually be installed on Linux and HPC systems using the appropriate package manager command (\ ``yum``\ , ``dnf``\ , or something similar). If you do not have administrative access to your machine, you may require the assistance of a system administrator. Once this installation has been performed, make sure to note any instructions for initializing your environment to use HDF5 (local setup scripts, ``module`` commands, etc.).

If you prefer to build the HDF5 library from source code (for example, to use a more recent release of the software), instructions may be found :doc:`here <build_hdf5>`.

----

Python and additional modules
-----------------------------

The ``kaiju`` software was developed assuming the use of **Python 3.8 or later**.

The installation of additional Python modules is described :doc:`here <install_python>`.

----

Git LFS (Large File Storage)
----------------------------

Our repository uses ``git lfs`` to efficiently store large binary files, so you must have ``git lfs`` installed and initialized before you clone our repository.

Many systems already have ``git lfs`` installed as part of their ``git`` installation. If this is the case you only need to run a one-time command so that ``git lfs`` can configure your ``git`` installation:

.. code-block::

   #!shell
   git lfs install


If your system does not have ``git lfs`` already installed, here are some example instructions for how to install it, ending with the same command as above: `Install git lfs <https://github.com/git-lfs/git-lfs/wiki/Installation>`_. You may need the assistance of a system administrator to perform this installation.

If you set up ``git lfs`` after already cloning the repository, you may retrieve the ``lfs`` files with

.. code-block::

   #!shell
   git lfs fetch
   git lfs pull


Note: To use ``git lfs`` on NASA HEC systems such as ``pleiades``\ , you need to load the ``pkgsrc`` module with ``module load pkgsrc``.

----

Cloning the ``kaiju`` repository
------------------------------------

NOTE: These instructions assume you have been given access to the BitBucket repository for the ``kaiju`` software. It further assumes that you have installed an SSH keypair for BitBucket access. If you need help with either of these steps, please see a CGS team member.

The ``kaiju`` repository can be cloned into any convenient location. In the following example, the code is cloned into the home directory.

.. code-block::

   #!shell
   cd $HOME
   git clone git@bitbucket.org:aplkaiju/kaiju.git


----

``cmake``
-------------

The ``kaiju`` software is built using the ``cmake`` system. This software is usually preinstalled on Linux and HPC systems. If you need to install it, instructions can be found `here <https://cmake.org/>`_.

----

Next steps
----------

Once these prerequisites are in place, you can proceed to `building the ``kaiju`` code <../building/build.md>`_.
