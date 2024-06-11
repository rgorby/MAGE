
Building a python environment for ``kaiju`` on ``pleiades``
===================================================================

Introduction
------------

This file describes how to set up a python 3.8-based ``python`` environment on ``pleiades`` which can be used to run the ``kaiju`` code.

These steps assume that ``pip`` (only) is used for installing additional packages, and packages are only taken from PyPI - no other repositories are used.

These instructions show how to install the ``miniconda`` python distribution into your home directory.

Building the python environment
-------------------------------

.. code-block::

   #!shell
   # Download the installer.
   cd $HOME
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

   # Run the installer.
   # Install into $HOME/pleiades/miniconda3, use all defaults.
   bash ./Miniconda3-latest-Linux-x86_64.sh

   # Make sure the conda configuration is loaded.
   source $HOME/.bashrc

   # Turn off auto-activation of base environment.
   conda config --set auto_activate_base false

   # Update everything to latest version.
   conda update --all

   # Now create the environment for kaiju, specifying only the python version.
   conda create -n kaiju-3.8 python=3.8

   # Activate the new environment.
   conda activate kaiju-3.8

   # Update everything in this environment to the latest version.
   conda update --all

   # Install packages required by kaipy (see below for contents of requirements.txt).
   pip install -r requirements.txt


The file ``requirements.txt`` is a text file containing:

.. code-block::

   #!shell
   ai.cdas
   alive_progress
   # cartopy
   cdasws
   cdflib
   configparser
   dataclasses_json
   h5py
   jinja2
   matplotlib
   netCDF4
   progressbar
   pyhdf
   pytest
   scipy
   spacepy
   sunpy
   xarray


If you have already cloned the ``kaiju`` repository, you will find a copy of this file as:

.. code-block::

   $KAIJUHOME/kaipy/requirements.txt


You can manually comment out (using the ``#`` syntax) the line for ``cartopy`` on ``pleiades``. If you do not comment out ``cartopy``\ , you will get errors during the installation process.

Using the python environment
----------------------------

Once your python environment is created, you must activate it for use with the ``kaiju`` software:

.. code-block::

   #!shell
   conda activate kaiju-3.8   # or whatever name you used for the python environment
