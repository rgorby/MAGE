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

For this example, we will use run the code on ``pleiades``, and use the default ``BASIC`` mode, which requires the minimum amount of input from the user. At each prompt, you can either type in a value, or hit the ``Return`` key to accept the default value (shown in square brackets at the end of the prompt). To get started, run ``makeitso.py`` with no arguments:

.. code-block:: bash

    $ source ~/local/cdf/3.9.0/bin/definitions.B
    $ makeitso.py
    Name to use for PBS job(s) [geospace]:
    Do you have an existing boundary condition file to use? (Y|N) [N]:
    Start date for simulation (yyyy-mm-ddThh:mm:ss) [2016-08-09T09:00:00]:
    Stop date for simulation (yyyy-mm-ddThh:mm:ss) [2016-08-09T11:00:00]:
    Do you want to split your job into multiple segments? (Y|N) [N]:
    GAMERA grid type (D|Q|O|H) [Q]:
    Name of HPC system (derecho|pleiades) [pleiades]:
    PBS account name [ewinter]:
    Run directory [.]:
    Path to kaiju installation [/home3/ewinter/cgs/aplkaiju/kaiju-private/development/kaiju-private]:
    Path to kaiju build directory [/home3/ewinter/cgs/aplkaiju/kaiju-private/development/kaiju-private/build_mpi]:
    PBS queue name (low|normal|long|debug|devel) [normal]:
    WARNING: You are responsible for ensuring that the wall time is sufficient to run a segment of your simulation!
    Requested wall time for each PBS job segment (HH:MM:SS) [01:00:00]:
    (GAMERA) Relative path to HDF5 file containing solar wind boundary conditions [bcwind.h5]:
    (VOLTRON) File output cadence in simulated seconds [60.0]:
    After these inputs, the script fetches data from CDAWeb for the specified time range to use in the solar wind boundary condition file, and generates XML and PBS files for the run, as well as a grid file for use in the model.

You should see output similar to this:

.. code-block:: bash

    Generating Quad LFM-style grid ...

    Output: lfmQ.h5
    Size: (96,96,128)
    Inner Radius: 2.000000
    Sunward Outer Radius: 30.000000
    Tail Outer Radius: 322.511578
    Low-lat BC: 45.000000
    Ring params:
    <ring gid="lfm" doRing="T" Nr="8" Nc1="8" Nc2="16" Nc3="32" Nc4="32" Nc5="64" Nc6="64" Nc7="64" Nc8="64"/>

    Writing to lfmQ.h5
    /home3/ewinter/miniconda3/envs/kaiju-3.8/lib/python3.8/site-packages/spacepy/time.py:2367: UserWarning: Leapseconds may be out of date. Use spacepy.toolbox.update(leapsecs=True)
    warnings.warn('Leapseconds may be out of date.'
    Retrieving f10.7 data from CDAWeb
    Retrieving solar wind data from CDAWeb
        Using Bx fields
    Bx Fit Coefficients are  [-3.78792744 -0.77915822 -1.0774984 ]
    Saving "OMNI_HRO_1MIN.txt_bxFit.png"
    Converting to Gamera solar wind file
        Found 21 variables and 120 lines
        Offsetting from LFM start ( 0.00 min) to Gamera start ( 0.00 min)
    Saving "OMNI_HRO_1MIN.txt.png"
    Writing Gamera solar wind to bcwind.h5
    Reading /home3/ewinter/cgs/aplkaiju/kaipy-private/ewinter-supermag_updates/kaipy-private/kaipy/rcm/dktable
    Reading /home3/ewinter/cgs/aplkaiju/kaipy-private/ewinter-supermag_updates/kaipy-private/kaipy/rcm/wmutils/chorus_polynomial.txt
    Dimension of parameters in Chorus wave model, Kp: 6 MLT: 97 L: 41 Ek: 155
    Wrote RCM configuration to rcmconfig.h5


    Template creation complete!


    The PBS scripts ['./geospace-00.pbs'] have been created, each with a corresponding XML file. To submit the jobs with the proper dependency (to ensure each segment runs in order), please run the script geospace_pbs.sh like this:
    bash geospace_pbs.sh

You should see the following new files in your run directory:

.. code-block:: bash

    $ ls
    bcwind.h5        geospace.json    OMNI_HRO_1MIN.txt_bxFit.png
    geospace-00.pbs  geospace_pbs.sh  OMNI_HRO_1MIN.txt.png
    geospace-00.xml  lfmQ.h5          rcmconfig.h5

The image files are summaries of the CDAWeb data used in the initial condition file (``bcwind.h5``). Those plots should look like:

.. image:: Bx_fit.png

.. image:: sw.png

Finally, submit the model run using the script generated by ``makeitso.py``. You will see the resulting PBS job ID.

.. code-block:: bash

    $ bash geospace_pbs.sh
    21602584.pbspl1.nas.nasa.gov

Once the job is started in the queue, it should take about XX to run. When complete, you will see the following in your run directory:

.. code-block:: bash

    $ ls
    XXX

Now perform a quick visualization of the results from your model using the ``msphpic.py`` script, provided in the ``kaipy`` package.

.. code-block:: bash

    $ msphpic.py -id geospace

This script will create a file called ``XXX.png``, which should look like this:

.. image:: arg1