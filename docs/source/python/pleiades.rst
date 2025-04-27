Creating a Python environment for ``kaiju`` on ``pleiades``
===========================================================


Introduction
------------

This file describes how to set up a python ``conda`` environment on
``pleiades`` which can be used to run the ``kaiju`` code.

These steps assume that ``pip`` (only) is used for installing additional
packages, and packages are only taken from PyPI.

These instructions assume you are using the ``bash`` shell.

.. note:: These instructions are designed for novice users. If you are
    comfortable with building and managing ``conda``-based environments, feel
    free to build your own environment using the NAS-provided
    `conda <https://www.nas.nasa.gov/hecc/support/kb/managing-and-installing-python-packages-in-conda-environments_627.html>`_
    software.


Building the python environment
-------------------------------

To create a Python environment for ``kaiju``, we first install Python, then
create the environment, then populate it with the required Python packages.

.. code-block:: bash


    # Download the Miniconda installer.
    cd $HOME
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86.sh

    # Run the installer.
    # Install into $HOME/miniconda3, use all defaults.
    bash ./Miniconda3-latest-Linux-x86.sh

    # Make sure the shell is properly configured.
    # Replace with rc file for your shell.
    source $HOME/.bashrc

    # Update everything to latest version.
    conda activate base
    conda update --all

    # Now create the environment for kaiju, specifying only the python
    # version.
    conda create -n kaiju-3.8 python=3.8

    # Activate the new environment.
    conda activate kaiju-3.8

    # Install the kaipy software.
    # If you are using kaipy via pip:
    pip install kaipy
    # OR
    # If you are using a clone of the kaipy repository:
    pip install -r /path/to/kaipy/requirements.txt


Using the python environment
----------------------------

Once your python environment is created, you must activate it for use with the
``kaiju`` software, and then enable the ``kaipy`` software:

.. code-block:: bash

    conda activate kaiju-3.8
    source /path/to/kaipy/kaipy/scripts/setupEnvironment.sh
    # or setupEnvironment.csh if using tcsh or a similar shell.
