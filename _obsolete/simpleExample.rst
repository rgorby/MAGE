A simple example
================

The instructions below walk you through the process of running a simple magnetosphere model to test your build of the ``kaiju`` code.

Preparing to run a ``kaiju`` model
----------------------------------

To set up your environment to run the ``kaiju`` software, the following steps are required.

    1. Load the same modules that you loaded when you built the ``kaiju`` software. For example, on `derecho`, you would run the following commands:

    .. code-block:: bash

        $ module --force purge
        $ module load ncarenv/23.06
        $ module load cmake/3.26.3
        $ module load craype/2.7.20
        $ module load intel/2023.0.0
        $ module load ncarcompilers/1.0.0
        $ module load cray-mpich/8.1.25
        $ module load hdf5-mpi/1.12.2

    2. *Source* (not *run*) the environment setup script for the ``kaiju`` software. For example, if the root of your ``kaiju`` repository clone is at ``$HOME/kaiju``, then you would run:

    .. code-block:: bash

        $ source $HOME/kaiju/scripts/setupEnvironment.sh

Running a simple magnetosphere problem
--------------------------------------

The ``kaiju`` software needs several files in order to run. The detailed steps for creating these files have been combined into a script called ``makeitso.py``. The script is provided in the ``kaiju`` code repository. You can see the options supported my ``makeitso.py`` by running it with the ``--help`` or ``-h`` command-line option.

.. code-block:: bash

    $ makeitso.py --help
    usage: makeitso.py [-h] [--clobber] [--debug] [--mode MODE] [--options_path OPTIONS_PATH] [--verbose]

    Interactive script to prepare a MAGE magnetosphere model run.

    optional arguments:
      -h, --help            show this help message and exit
      --clobber             Overwrite existing options file (default: False).
      --debug, -d           Print debugging output (default: False).
      --mode MODE           User mode (BASIC|INTERMEDIATE|EXPERT) (default: BASIC).
      --options_path OPTIONS_PATH, -o OPTIONS_PATH
                            Path to JSON file of options (default: None)
      --verbose, -v         Print verbose output (default: False).