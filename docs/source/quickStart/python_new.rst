
Python modules
==============

Introduction
------------

The `standard Python distribution <https://www.python.org/>`_ includes the standard Python library. However, several additional Python packages are required for running the ``kaiju`` software. These additional modules can be installed in one of 3 standard ways (there are others, but these are by far the most common):


* Build from a source tarball.
* Use the ``pip`` package manager.
* Use the ``conda`` package manager.

Code built from a source tarball provides the most flexibility in your installation, but also requires the most time, as well as detailed knowledge of build procedures. We will only use this method when the required modules are not available via ``pip`` or ``conda``.

The ``pip`` and ``conda`` tools provide the same function - Python package management. The major difference (other than command syntax) is the default software source used by the tools. For ``pip``\ , packages are (by default) installed from the `Python Package Index <https://www.pypi.org>`_\ , which is the standard package repository for the Python community. The ``conda`` tool by default installs packages from repositories managed by `\ ``anaconda.org`` <https://anaconda.org/>`_\ , managed by `Anaconda Inc. <https://www.anaconda.com/>`_ as part of its Anaconda Python distribution. The ``conda`` package manager is also available using the minimal distribution provided by `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_. We strongly prefer to use the miniconda-based installation on all platforms, since it gives the user complete control over the precise Python package versions in use. However, most other python environments should "just work". For example, users have reported success using this procedure with the python "module" system provided on pleiades.

These instructions will assume use of the ``miniconda`` distribution of python. Package management will use ``pip``.

Installing ``miniconda``
----------------------------

This step may be skipped if you already have a working ``miniconda`` or Anaconda installation.


#. 
   Download the Miniconda installer for your computer from the `Miniconda web page <https://docs.conda.io/en/latest/miniconda.html>`_.

#. 
   Run the installer:

    #!shell
    # Installer name will be different for non-Mac systems.
    bash ./Miniconda3-latest-MacOSX-x86_64.sh

Accept defaults for all prompts from the installer.


#. 
   Source your ``~/.bashrc`` to activate your installation:

    #!shell
    source ~/.bashrc

Creating a virtual environment for the ``kaiju`` software
-------------------------------------------------------------

We strongly suggest creating a separate python virtual environment to use with the ``kaiju`` software. These instructions assume use of python 3.8. Later versions of python should also work.


#. 
   Create a new virtual python environment using python 3.8:

    #!shell
    conda create -n kaiju-3.8 python=3.8

You can replace ``kaiju-3.8`` with any name you wish to use for your environment. You may replace ``3.8`` with your desired python version.


#. 
   Activate the new environment:

    #!shell
    conda activate kaiju-3.8

#. 
   Update all software in the environment:

    #!shell
    conda update --all

Installing additional Python packages
-------------------------------------

The Python part of the ``kaiju`` software requires many other python packages to function. Once you have activated your virtual environment for ``kaiju``\ , you can install the required packages as follow:

.. code-block::

   #!shell
   pip install ai.cdas
   pip install alive_progress
   pip install cartopy
   pip install cdasws
   pip install cdflib
   pip install configparser
   pip install dataclasses_json
   pip install h5py
   pip install matplotlib
   pip install netCDF4
   pip install progressbar
   pip install pyhdf
   pip install pytest
   pip install scipy
   pip install spacepy
   pip install sunpy
   pip install xarray


Alternatively, you can install these packages with a single command using a package list provided as part of the ``kaiju`` software:

.. code-block::

   #!shell
   pip install -r $KAIJU_HOME/kaipy/requirements.txt


These commands install the latest versions of the requested modules, which should be fine for the latest version of the ``kaiju`` software.
