
Building a python environment for ``kaiju`` on ``derecho``
==================================================================

Introduction
------------

This file describes how to set up a python 3.8-based ``conda`` environment on ``derecho`` which can be used to run the ``kaiju`` code.

NOTE: Since the same home directory is used for ``casper``\ , ``cheyenne``\ , and ``derecho``\ , the ``conda`` installations for each must go in separate subdirectories. However, unless otherwise directed, all of the installations will share the same ``pip`` cache, ``~/.condarc``\ , etc.. That means if you build a wheel for a compiled module on one system, that wheel will be found and used on all of the other systems. This can cause problems since the CPU architectures and default module environments on the 3 systems differ.

These steps assume that ``pip`` (only) is used for installing additional packages, and packages are only taken from PyPI (and the default ``conda`` repository, if needed) - no other repositories are used.

.. code-block::

   #!shell
   # Download the installer.
   cd $HOME/derecho
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

   # Run the installer.
   # Install into $HOME/derecho/miniconda3, use all defaults.
   # Run conda init, then update ~/.bashrc to only run this on derecho.
   bash ./Miniconda3-latest-Linux-x86_64.sh

   # NOTE: This installation *should* run "conda init", which *should*
   # add the conda setup code to ~/.bashrc.

   # NOTE: This installation creates ~/.conda/environments.txt, but
   # nothing else. This file just contains the path to the miniconda3
   # installation directory.

   # Turn off auto-activation of base environment.
   # NOTE: This command creates ~/.condarc.
   conda config --set auto_activate_base false

   # Update everything to latest version.
   conda update --all

   # NOTE: This creates the directory ~/.cache/conda.

   # Now create the environment for kaiju, specifying only the python version.
   conda create -n kaiju-3.8 python=3.8

   # NOTE: This adds the path to the new environment in
   # ~/.conda/environments.txt.

   # Activate the new environment.
   conda activate kaiju-3.8

   # IMPORTANT: Use a pip cache specific to this machine.
   export PIP_CACHE_DIR=$HOME/derecho/pip_cache

   # Install the latest versions of required packages.
   conda update --all

   # This module is needed to install cartopy.
   module load geos/3.9.1

   # Install packages required by kaipy.
   pip install -r requirements.txt

   # Save the environment description.
   conda list >& kaiju-3.8-derecho-YYYYMMDD.env


Where ``requirements.txt`` is a text file containing:

.. code-block::

   #!shell
   ai.cdas
   alive_progress
   cartopy
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


Using the python environment
----------------------------

Once your python environment is created, you must activate it for use with the ``kaiju`` software:

.. code-block::

   #!shell
   conda activate kaiju-3.8   # or whatever name you used for the python environment


If you want environment to be available for Jupyter interface you need to install ipykernel into your evironment

.. code-block::

   conda activate kaiju-3.8
   mamba install ipykernel
