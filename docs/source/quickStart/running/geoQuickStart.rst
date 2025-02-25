Geospace Quick Start Guide
==========================

.. note::

    These instructions assume you are running on ``derecho``. The instructions
    for ``pleiades`` differ only in the list of loaded modules.

Preparing to run a ``kaiju`` model
----------------------------------

To set up your environment to run the ``kaiju`` software, the following steps
are required.

    1. Load the same modules that you loaded when you built the ``kaiju``
    software. For example, on ``derecho``, you would run the following commands:

    .. code-block:: bash

        $ module --force purge
        $ module load ncarenv/23.06
        $ module load cmake/3.26.3
        $ module load craype/2.7.20
        $ module load intel/2023.0.0
        $ module load ncarcompilers/1.0.0
        $ module load cray-mpich/8.1.25
        $ module load hdf5-mpi/1.12.2

    2. *Source* (not *run*) the environment setup scripts for the ``kaiju``
    software *and* the ``kaipy`` software. For example, if the root of your
    ``kaiju`` repository clone is at ``$HOME/kaiju``, and the root of your
    ``kaipy`` repository clone is at ``$HOME/kaipy``, then you would run:

    .. code-block:: bash

        $ source $HOME/kaiju/scripts/setupEnvironment.sh
        $ source $HOME/kaipy/kaipy/scripts/setupEnvironment.sh

    These scripts set several environment variables, including ``KAIJUHOME``,
    which refers to the root directory of your ``kaiju`` software, and
    ``KAIPYHOME``, which refers to the root directory of your ``kaipy``
    software.

Create the grid file
--------------------

The ``kaiju`` software requires an HDF5 file which provides the coordinates of
the vertices of the grid over which the model is to be run. These files can be
generated easily using the ``genLFM.py`` script, which is distributed with
``kaipy``. The script can generate LFM (Lyon-Fedder-Mobarry) grids of the
magnetosphere in one of four resolutions: (D)ouble, (Q)uad, (O)ct, (H)ex
(listed in order of increasing grid resolution). For example, the following
command will generate a Q-resolution LFM grid file:

.. code-block:: bash

    $ genLFM.py -gid Q
    Generating Quad LFM-style grid ...

    Output: lfmQ.h5
    Size: (96,96,128)
    Inner Radius: 2.000000
    Sunward Outer Radius: 28.000106
    Tail Outer Radius: 301.011944
    Low-lat BC: 45.000000
    Ring params:
    <ring gid="lfm" doRing="T" Nr="8" Nc1="8" Nc2="16" Nc3="32" Nc4="32" Nc5="64" Nc6="64" Nc7="64" Nc8="64"/>

    Writing to lfmQ.h5

Create the solar wind boundary condition file
---------------------------------------------

The simplest way to generate a solar wind boundary condition file is to use
the ``cda2wind.py`` script, which is part of ``kaipy``. The script downloads
data from `CDAWeb <https://cdaweb.gsfc.nasa.gov/>`_ and converts it into the
format required by the ``kaiju`` software. For example, to request solar wind
data for the 2-hour period beginning at 0900 UT on 9 August 2016, you would
enter the command:

.. code-block:: bash

    $ cda2wind.py -t0 2016-08-09T09:00:00 -t1 2016-08-09T11:00:00

This command will create a file called ``bcwind.h5`` in the current directory
which contains the solar wind data needed for the inner boundary conditions of
the model.

Create the XML input file
-------------------------

The input deck is an XML file that specifies various parameters to be read at
runtime. A typical XML file might look like this:

.. code-block:: xml

    <?xml version="1.0"?>
    <Kaiju>
        <Gamera>
            <sim H5Grid="lfmQ.h5" doH5g="T" icType="user" pdmb="0.75" rmeth="8C" runid="geospace" />
            <floors dFloor="1.0e-4" pFloor="1.0e-6" />
            <timestep doCPR="T" limCPR="0.20" />
            <restart doRes="F" nRes="-1" resID="geospace" />
            <physics Ca="10.0" doBoris="T" doMHD="T" />
            <ring doRing="T" gid="lfm" />
            <wind tsfile="bcwind.h5" />
            <source doBounceDT="T" doSource="T" doWolfLim="T" nBounce="1.0" />
            <iPdir N="4" bcPeriodic="F" />
            <jPdir N="4" bcPeriodic="F" />
            <kPdir N="1" bcPeriodic="T" />
            <coupling blockHalo="T" />
        </Gamera>
        <REMIX>
            <conductance apply_cap="T" const_sigma="F" doStarlight="T" />
            <precipitation aurora_model_type="LINMRG" beta="0.2" doAuroralSmooth="F" />
        </REMIX>
        <VOLTRON>
            <time tFin="7200.0" />
            <spinup doSpin="T" tSpin="7200.0" />
            <output dtOut="60.0" tsOut="300.0" />
            <coupling doAsyncCoupling="F" doQkSquish="T" dtCouple="5.0" imType="RCM" qkSquishStride="2" />
            <restart dtRes="1800.0" />
            <imag doInit="T" />
            <helpers doSquishHelp="T" numHelpers="2" useHelpers="T" />
        </VOLTRON>
        <CHIMP>
            <units uid="EARTHCODE" />
            <fields grType="lfm" />
            <domain dtype="MAGE" />
            <tracer epsds="0.05" />
        </CHIMP>
        <RCM>
            <rcmdomain domType="ELLIPSE" />
            <ellipse isDynamic="T" xSun="12.5" xTail="-15.0" yDD="15.0" />
            <grid HiLat="75.0" LowLat="30.0" />
            <plasmasphere doRefill="T" initKp="5" isDynamic="T" tAvg="60.0" />
        </RCM>
    </Kaiju>

Some of the important parameters are:

``tFin="7200.0"``
    The run will proceed for 7200 simulated seconds (2 simulated hours).

``dtOut="60.0"``
    Screen output will be produced every 60 simulated seconds (1 simulated
    minute

``dtRes="1800.0"``
    Create restart files every 1800 simulated seconds (30 simulated minutes).

``dtCouple="5.0"``
    Perform REMIX-GAMERA coupling every 5 simulated seconds (). Note that
    ``dtCouple`` is a required input for ionospheric coupling.

``runid="geospace"``
    The run ID of the model resuslts.

``H5Grid="lfmQ.h5"``
    Name of the grid file created with ``genLFM.py``.

``doRes="F"``
    This job will not begin with a restart file.

``<ring>``
    ``gid`` tells the model which type of grid is being used. It is supposed
    to be consistent with the input grid file. Ring average technique is
    implemented if ``doRing`` is ``T``. The number of parameters for ring
    average is set in ``Nr``, and each parameter is listed as ``Nc1``,
    ``Nc2``, ... Usually there are four parameters for a double resolution
    grid and 8 parameters for a quad resolution grid.

``tsfile="bcwind.h5"``
    Specify the name of the solar wind file to be used.

.. important:: 
    The initial state (a pure dipole) with the conductance model turned on can
    sometimes be erratic.  It's better to spin up for a while using constant
    conductance then do a restart where the conductance model is turned on.

Create the PBS script
---------------------

Executable and run script
~~~~~~~~~~~~~~~~~~~~~~~~~

Compile the target ``voltron.x`` and then just run it with a PBS script.  You
should probably request a full node and run with 72 threads.  Below is the
PBS script that I use, but you might need to make changes to it for your
setup.  This is assuming an XML file called "cmiD.xml", that you have a
preset module savelist, and that voltron.x is in your path.  Probably
some other stuff too.  I ran it with

.. code-block:: bash

   qsub -v CHIMPEXE="voltron.x" -N cmiD RunCHIMP.pbs

and it worked.

.. code-block:: bash

   #!bash
   #!/bin/bash
   #PBS -A P28100045
   #PBS -N KAIJU
   #PBS -j oe
   #PBS -q regular
   #PBS -l walltime=12:00:00
   #PBS -l select=1:ncpus=72:ompthreads=72

   #Example usage
   #qsub -v CHIMPEXE="slice.x" -N cmiD RunCHIMP.pbs
   #Module savelist "kaiju"
   #Currently Loaded Modules:
   #  1) git/2.9.5    (H)   4) impi/2018.4.274       7) python/2.7.16
   #  2) intel/18.0.5       5) ncarenv/1.3           8) cmake/3.14.4
   #  3) hdf5/1.10.5        6) ncarcompilers/0.5.0

   export EXE=${CHIMPEXE:-"slice.x"}
   export RUNID=${PBS_JOBNAME}

   source ~/.bashrc
   module restore kaiju

   module list
   hostname
   date
   export OMP_NUM_THREADS=72
   export KMP_STACKSIZE=128M
   export JNUM=${PBS_ARRAY_INDEX:-0}
   echo "Running $EXE"
   ${EXE} ${RUNID}.xml ${JNUM} > ${RUNID}.${JNUM}.out
   date

Running the model
-----------------

Visualization
-------------

Data can be read into VisIt or Paraview using the ``kaiju/scripts/genXDMF.py``
script which will create an XDMF file which either of those viewers can
natively read.

Alternatively, you can try ``kaiju/scripts/msphpic.py`` script which uses the
Python post-processing routines to generate a multi-panel plot of a
magnetosphere run.  It uses Python libraries that you may or may not have,
I don't know your life.

Creating gamera videos in parallel via slurm job submission and gamsphVid.py.

.. code-block:: bash

    #!/bin/bash -l
    #SBATCH -J VidMHD
    #SBATCH --output=%x_%A_%a.out
    #SBATCH -t 24:00:00
    #SBATCH -A UJHB0010
    #SBATCH -p dav
    #SBATCH --array=1-36
    #Defaults
    export BASE="/glade/u/home/skareem/Work/gamercm/Data/"
    export RUNID="vapD9"

    export TS="60"
    export TE="1480"
    export DT="60"

    export ARG=""
    export ODIR="VidMHD"
    export JNUM=${SLURM_ARRAY_TASK_COUNT:-"1"}
    export JID=${SLURM_ARRAY_TASK_ID:-"1"}

    echo "My JobID: " $SLURM_ARRAY_TASK_ID " of " $SLURM_ARRAY_TASK_COUNT
    export IDIR="${BASE}${RUNID}"
    source ~/.bashrc
    setdav
    gamsphVid.py -d $IDIR -ts $TS -te $TE -dt $DT -Nblk $JNUM -nID $JID -o $ODIR $ARG

Modify and launch the above script using ``sbatch ABOVEFILE.s``:

.. code-block:: bash

    #sbatch array#

Specifies how many parallel jobs you want (36 here)

.. code-block:: bash

    export TS="60"
    export TE="1480"
    export DT="160"

are the time domain in which you want to create movies.  TS is the start time
(in minutes) of your run starting from T=0.  TE is time end (in minutes).  DT
is the time (in seconds) between slices.  DT will attempt to find the cloest
timestep, so if outputs are 60s but you input DT=15s, then you get duplicated
frames.

.. code-block:: bash

    source ~/.bashrc
    setdav

Set up your environment variables.  These are Kareem specific, so modify to
match whatever your Casper environment is.

MPI Magnetosphere
-----------------

Running a basic Gamera case with MPI support requires three things:

#. Building the MPI version of the gamera executable
#. Modifying case XML to supply additional MPI decomposition information
#. Modifying the submission script to request multiple nodes and use mpirun

Modifying Case XML
~~~~~~~~~~~~~~~~~~

Three additional lines are required in the case XML file when running with
MPI decomposition. And these lines are ignored by the non-MPI version of
Gamera, so you can safely leave them in the XML (if you want to) when not
using MPI.

In the ``<Gamera>```\ "`` section of the XML, one line is required for each
dimension that tells how many regions that dimension is decomposed into, and
whether that dimension is periodic. Here is an example where the case is
decomposed into 4 regions each along the I and J axes, and 2 regions along the
K axis. The I and J axes are not periodic, but the K axis is.

.. code-block:: xml

    <Gamera>
    ...
        <iPdir N="4" bcPeriodic="F"/>
        <jPdir N="4" bcPeriodic="F"/>
      < kPdir N="2" bcPeriodic="T"/>
    ...
    </Gamera>

(Corrections: Only k=1 works with the magnetosphere because the ringav code
needs all the k direction data on one rank. Change to kPdir N="1". Accordingly
the PBS job script should use:

.. code-block:: bash

    #PBS -l select=8:ncpus=36:mpiprocs=2:ompthreads=18+1:ncpus=36:mpiprocs=1:ompthreads=36

with the part up to the plus sign allocating for gamera proper, following the
plus sign allocates for voltron.

Modifying Job Submission Script
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Assuming that you have a job submission script suitable for running serial
jobs, here is how to modify it to run MPI gamera. These examples will be
designed to work with Cheyenne, but can be adapted to most clusters.

First, you need to request an appropriate number of nodes and tell the job
submission system how many MPI ranks should be created per node. Continuing
the example above, this case is decomposed 4 times along the I and J
dimensions, and twice along the K dimension. So this case will need a total of
4*4*2=32 MPI ranks. For this example we will assign one MPI rank to each
physical processor/socket, which provides a reasonable balance between
performance and cost. Each of Cheyenne's compute nodes has two
processors/sockets, so this means that each compute node will receive two MPI
ranks. The original, serial, resource request line from the job submission
script looked like this:

.. code-block:: bash

    #PBS -l select=1:ncpus=72:ompthreads=72

That line requests all 72 cpus on a single compute node, which is perfect for
a single process. We want to create a total of 32 processes (one per MPI rank)
spread across 16 compute nodes (2 MPI ranks per node), with each MPI rank
getting half of the compute node's compute resources. That looks like this:

.. code-block:: bash

    #PBS -l select=16:ncpus=36:mpiprocs=2:ompthreads=36

The only other line we need to change is the one that calls the executable
and starts the simulation. In the serial case that looked like this:

.. code-block:: bash

    ${EXE} ${RUNID}.xml ${JNUM} > ${RUNID}.${JNUM}.out

That command literally calls the executable and passes it the input XML file.
Instead, we now need to use a helper application called mpirun, which will
call our executable for us:

.. code-block:: bash

    mpirun ${EXE} ${RUNID}.xml ${JNUM} > ${RUNID}.${JNUM}.out

Using MPT
---------

If you are running with the MPT mpi library, the submission script will
require some additional modifications, described in a dedicated wiki page
:doc:`Running with MPT </userGuide/magnetosphere/runningWithMPT>`.

Generating solar wind file manually
-----------------------------------

Solar wind data can be defined using an HDF-5 file.  This is specified in the
XML input deck under "Gamera/wind/tsfile" as seen in the example XML file
above.

Solar wind HDF5 file input units

* Time [s]
* Density [#/cc]
* Velocity [m/s]
* Pressure [nPa]
* Magnetic field [nT]

Gamera will read the HDF file and convert to its internal code units (see
normalization above).  Below is an example of a Python script to generate a
solar wind file.

The repository includes a Jupyter Notebook,
``analysis/notebooks/Tutorial/Kaiju-IdealSW.ipynb`` that provides an
interactive way to create your own solar wind file.  You can also modify the
python code below to create a solar wind file.

.. code-block:: python

    #!python

    import h5py
    import numpy as np
    import sys
    import kaiTools as kT

    TW = 1.0e+5     #Default temperature, K
    nW = 5          #Default density, #/cc
    VxW = 400.0     #Default wind, km/s
    f107val = 100.0 #Default f10.7 flux
    tilt = 0.0      #Default dipole tilt, radians
    mjd0 = 58767.0  #Default MJD, set for 2019-10-11 00:00:00

    Bx0 = 0.0 #Default Bx offset for planar front, keep at zero
    ByC = 0.0 #Default By coefficient used to calculate Bx, include if want tilted field
    BzC = 0.0 #Default Bz coefficient used to calculate Bx, include if want tilted field

    fOut = "bcwind.h5"

    #Time bounds [hours]
    tMin = 0.0
    tMax = 6.0
    dt = 60.0 #Cadence [s]

    SimT = (tMax-tMin)*60.0*60.0
    NumT = np.int( np.ceil(SimT/dt)+1 )

    print("Generating %d slices, T=[%5.2f,%5.2f]"%(NumT,tMin,tMax))

    T = np.linspace(tMin,tMax,NumT)
    D = np.zeros(NumT)
    Temp = np.zeros(NumT)
    Vx = np.zeros(NumT)
    Vy = np.zeros(NumT)
    Vz = np.zeros(NumT)
    Bx = np.zeros(NumT)
    By = np.zeros(NumT)
    Bz = np.zeros(NumT)
    f107 = np.zeros(NumT)
    ThT = np.zeros(NumT)
    mjd = np.zeros(NumT)
    symh = np.zeros(NumT)

    tWin = 1.0 #Window times [hr]
    for i in range(NumT):
        t = T[i] #Time in hours
        if (t <= tWin):
            D[i] = nW
            Vx[i] = -VxW
            Temp[i] = TW
            f107[i] = f107val
            ThT[i] = tilt
            mjd[i] = mjd0 + T[i]/24.0
        elif (t <= 3*tWin):
            D[i] = nW
            Vx[i] = -VxW
            Temp[i] = TW
            Bz[i] = -5.0
            f107[i] = f107val
            ThT[i] = tilt
            mjd[i] = mjd0 + T[i]/24.0
        elif (t <= 6.0*tWin):
            D[i] = nW
            Vx[i] = -VxW
            Temp[i] = TW
            Bz[i] = +5.0
            f107[i] = f107val
            ThT[i] = tilt
            mjd[i] = mjd0 + T[i]/24.0
        else:
            D[i] = nW
            Vx[i] = -VxW
            Temp[i] = TW
            Bz[i] = -5.0
            f107[i] = f107val
            ThT[i] = tilt
            mjd[i] = mjd0 + T[i]/24.0

    #Write solar wind
    #t,D,V,Temp,B = [s],[#/cm3],[m/s],[K],[nT]

    oTScl = (60*60.0) #hr->s
    oDScl = 1.0
    oVScl = 1.0e+3 #km/s->m/s
    oTempScl = 1.0
    oBScl = 1.0

    #Compute Kp and Dst
    kp = kT.newellkp(time,Den,Vx,Vy,Vz,Bx,By,Bz)
    dst = kT.burtonDst(time,Den,Vx,Vy,Vz,Bx,By,Bz)

    with h5py.File(fOut,'w') as hf:
        hf.create_dataset("T" ,data=oTScl*T)
        hf.create_dataset("symh" ,data=dst)
        hf.create_dataset("Kp",data=kp)
        hf.create_dataset("D" ,data=oDScl*D)
        hf.create_dataset("Temp" ,data=oTempScl*Temp)
        hf.create_dataset("Vx",data=oVScl*Vx)
        hf.create_dataset("Vy",data=oVScl*Vy)
        hf.create_dataset("Vz",data=oVScl*Vz)
        hf.create_dataset("Bx",data=oBScl*Bx)
        hf.create_dataset("By",data=oBScl*By)
        hf.create_dataset("Bz",data=oBScl*Bz)
        hf.create_dataset("tilt",data=ThT)
        hf.create_dataset("f10.7",data=f107)
        hf.create_dataset("MJD",data=mjd)
        hf.create_dataset("Bx0",data=Bx0)
        hf.create_dataset("ByC",data=ByC)
        hf.create_dataset("BzC",data=BzC)
