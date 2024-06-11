
Building on ``derecho``
===========================

Introduction
------------

This page provides instructions for building the ``kaiju`` software on the ``derecho`` supercomputer.

Step 0: Prepare your environment
--------------------------------

Like most HPC systems, ``derecho`` uses the ``module`` system to manage the versions of software packages available to the user. A "module" is a collection of programs and libraries for a specific task, and "loading" the module adjusts the user environment (mostly by setting or updating environment variables) to make that module available to the user. When you log in to ``derecho``\ , the following modules are loaded by default:

.. code-block::

   #!shell
   module list

   Currently Loaded Modules:
     1) ncarenv/23.06 (S)   2) craype/2.7.20   3) intel/2023.0.0   4) ncarcompilers/1.0.0   5) cray-mpich/8.1.25   6) hdf5/1.12.2   7) netcdf/4.9.2


This default set of modules *cannot* be used to build the ``kaiju`` software. The correct module sets for the serial and MPI versions of ``kaiju`` will be listed in the instructions below.

Step 1: Build prerequisite libraries
------------------------------------

NASA CDF library
^^^^^^^^^^^^^^^^

The `NASA CDF (Common Data Format) library <https://cdf.gsfc.nasa.gov/>`_ is used in parts of the ``kaiju`` post-processing software when fetching spacecraft data from `CDAWeb <https://cdaweb.gsfc.nasa.gov/>`_. Prior to building the ``kaiju`` software, the CDF library must be built and installed, since it is not available as a module on ``derecho``.

Instructions for building and installing the CDF library on ``derecho`` are available `here <derecho_build_cdf>`_.

Step 2: Create a python environment
-----------------------------------

Most of the ``kaiju`` software for pre-processing, post-processing, and analysis is written in `Python <https://www.python.org/>`_. Python is available in many forms (or "distributions"), but we recommend use of the `Miniconda] distribution <https://docs.conda.io/en/latest/miniconda.html>`_ for simplicity and compactness.

Instructions for installing python and building a python environment on ``derecho`` are available `here <derecho_build_python>`_.

Step 3: Compile the ``kaiju`` software
------------------------------------------

The ``kaiju`` software can be built in serial or MPI versions. The serial version should be used when running the code on a single computer, such as a Mac laptop. The MPI version should be used when running on an HPC system (a supercomputer).

Instructions for building the serial and MPI versions of the ``kaiju`` software on ``derecho`` are available `here <derecho_build_kaiju>`_.
