.. role:: raw-html-m2r(raw)
   :format: html


Geospace Quick Start Guide
==========================

** Notes for Updating ** Want to replace this page the instructions on
steps to follow to do your own event simulation. Basic points, generate
solar wind file, decide on what components Gamera, Mix, RCM, T*GCM
you'll be using, decide on resolution, generate grid files, generate
batch script files, run the code, refer them to the post processing
steps.

**Note:** This quick start assumes you have completed the `build instructions <../building/build>`_.

**Note:** Throughout the descriptions ##$KAIJUHOME## refers to the
base directory of the `kaiju <https://bitbucket.org/aplkaiju/kaiju>`_
repository.

----

Magnetosphere Quick Start
=========================

This is a quick set of instructions to run a coupled Gamera-ReMIX (magnetosphere-ionosphere) run on Cheyenne.  It uses the executable "voltron.x"

[TOC]

OMP Magnetosphere
=================

Basic magnetosphere runs require five things


* Grid file: Collection of grid corners in HDF5 format
* Solar wind: Solar wind time series in HDF5 format
* Input deck: Run configuration parameters in XML format
* Executable: Duh (voltron.x)
* Run script: Job script to setup run for a single node on Cheyenne

Note, this is assuming the main repo has been downloaded to a base directory kaiju

Grid file
---------

The simplest way of generating is to use kaiju/scripts/genLFM.py to create a standard LFM-style grid of resolution D/Q/O/H for double/quad/oct/hex.  Note, you'll likely have to do this on Casper, the post-processing node.

.. code-block::

   [skareem@casper26:~/kaiju/scripts]$ genLFM.py -gid D
   Generating Double LFM-style grid ...

   Output: lfmD.h5
   Size: (48,48,64)
   Inner Radius: 2.000000
   Sunward Outer Radius: 28.000106
   Tail Outer Radius: 301.011944
   Low-lat BC: 45.000000
   Ring params:
   <ring gid="lfm" doRing="T" Nr="4" Nc1="8" Nc2="16" Nc3="32" Nc4="32"/>

   Writing to lfmD.h5

The grid generation spits parameters for ring average which will be added to the XML input deck.

Solar wind
----------

Gamera reads solar wind data using HDF5. The simplest way to generate an appropriate HDF5 file is to use the cda2wind.py script, which downloads data from the CDAS database between the specified timeframe and converts it into the correct format to be read in by Gamera. These scripts rely on python libraries that are listed in our `prerequisites <prerequisites>`_ so make sure you have them installed.

Solar Wind Generator
^^^^^^^^^^^^^^^^^^^^

The solar wind generator has numerous options and features that a described `here <solarWindGenerator>`_.

Converting LFM Solar Wind File
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An HDF5 file can also be created from an existing LFM-style solar wind input using the reWind2.py script. For an example see kaiju/examples/earthcmi/OMNI_HRO_1MIN_16772.txt_SW-SM-DAT. Note the GAMERA input of solar wind file requires more information than the density, velocity and IMF in LFM. The extra variables include F10.7, diple tilt, MJD (Julian date of each record), fitting coefficients of IMF (Bx0, ByC, BzC). Sound speed is not taken but the temperature is needed. 

The MJD is especially helpful when a different time step of solar wind input is needed. The default cadence of solar wind input is 1 min. In order to use higher rate solar wind input, we used to have to modify the init-fortran code in LFM. Here in the GAMERA world, the solar wind file MJD is the only variable to be modified for a flexible input cadence.

.. code-block::

   [skareem@casper26:~]$ reWind2.py OMNI_HRO_1MIN_16772.txt_SW-SM-DAT
   Reading LFM solar wind from OMNI_HRO_1MIN_16772.txt_SW-SM-DAT
       Found 11 variables and 1440 lines
       Offsetting from LFM start ( 0.00 min) to Gamera start ( 0.00 min)
   Writing Gamera solar wind to bcwind.h5

Input deck
----------

The input deck is an XML file that specifies various parameters to be read at runtime.  An example for a double resolution run is shown below.  Note that times used to be input in code units and the conversion is 63.8 seconds. Now we have changed the time variables in xml to be in seconds. This means

**NOTE: As of 2022-04-05, this XML does not conform to the current required XML format. It needs an enclosing ``<Kaiju>`` element, and perhaps other changes.**

.. code-block::

   #!xml
   <?xml version="1.0"?>
   <!-- Magnetosphere params, Voltron times in seconds -->
   <!-- MJD0 is modified Julian date of T=0 in solar wind file -->
   <VOLTRON>
     <time tFin="36000.0"/>
     <output dtOut="60.0" tsOut="100" doTimer="F"/>
     <restart dtRes="1800.0"/>
     <coupling dt="5.0"/>
   </VOLTRON>
   <Gamera>
     <sim runid="msphere" doH5g="T" H5Grid="lfmD.h5" icType="user" pdmb="1.0" pFloor="1.0e-8" dFloor="1.0e-6" rmeth="7UP"/>
     <restart doRes="F" resId="msphere" nRes="0"/>
     <physics doMHD="T" doBoris="T" Ca="10.0"/>
     <prob Rho0="0.2" P0="0.001"/>
     <ring gid="lfm" doRing="T" Nr="4" Nc1="8" Nc2="16" Nc3="32" Nc4="32"/>
     <wind tsfile="bcwind.h5"/>
   </Gamera>
   <!-- Remix params -->
   <REMIX>
     <grid Np="360" Nt="45" LowLatBoundary="45.0"/>
     <conductance pedmin="2.0" hallmin="1.0" sigma_ratio="3.0" const_sigma="True" ped0="5.0"/>
   </REMIX>

With these parameters the run will go for 10 hours (36000 s), outputting every 1 minute (dtOut=60.0 s) and writing restarts every 30 minutes (dtRes=1800.0 s). The ReMIX-Gamera coupling is done every 5 seconds (coupling dt=5.0). Note the coupling dt is a required input for the ionospheric coupling.

The runid gives the name of the output mhd/mix file. H5Grid is the name of the grid file that has to be present in the run directory. If the job to be submitted with this xml file is not restarting from a previously saved state, set doRes as False. Otherwise, to restart, set doRes as True and make sure the resFile is linked to the right restarting file. As of 23 April 2020 restarts are specified with ID/# instead of filename. Instead of restart/resFile, specify restart/resID and restart/nRes.
 The restart file msphere.Res.00005.h5 would be: 
    :raw-html-m2r:`<restart resId="msphere" nRes="5"/>` 
Specifying nRes="-1" will read the XXXXX symbolic link.

:raw-html-m2r:`<physics>` domain:

:raw-html-m2r:`<prob>` domain: 

:raw-html-m2r:`<ring>` domain: gid tells the model which type of grid is being used. It is supposed to be consistent with the input grid file. Ring average technique is implemented if doRing is True. The number of parameters for ring average is set in Nr, and each parameter is listed as Nc1, Nc2, ... Usually there are four parameters for double resolution grid and 8 parameters for quad resolution grid.

:raw-html-m2r:`<wind>` domain: tsfile takes the name of the solar wind file to be used.

Under the REMIX block, the Np and Nt gives the number of grid cells in Remix along longitude and latitude. The low latitude boundary is set at 45 deg latitude.

This setup uses constant conductance, an example of using the conductance model would replace the conductance block with this,

.. code-block::

   #!xml
   <conductance F107="100.0" pedmin="2.0" hallmin="1.0" sigma_ratio="3.0" const_sigma="False" ped0="10.0"/>

Note, that starting the initial state (a pure dipole) with the conductance model turned on can sometimes be erratic.  It's better to spin up for a while using constant conductance then do a restart where the conductance model is turned on.  But whatever, you do you.

Executable and run script
-------------------------

Compile the target "voltron.x" and then just run it with a PBS script.  You should probably request a full node and run with 72 threads.  Below is the PBS script that I use, but you might need to make changes to it for your setup.  This is assuming an XML file called "cmiD.xml", that you have a preset module savelist, and that voltron.x is in your path.  Probably some other stuff too.  I ran it with 

.. code-block::

   qsub -v CHIMPEXE="voltron.x" -N cmiD RunCHIMP.pbs

and it worked.

.. code-block::

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

Generating solar wind file manually
-----------------------------------

Solar wind data can be defined using an HDF-5 file.  This is specified in the XML input deck under "Gamera/wind/tsfile" as seen in the example XML file above.

 Solar wind HDF5 file input units


* Time [s]
* Density [#/cc]
* Velocity [m/s]
* Pressure [nPa]
* Magnetic field [nT]

Gamera will read the HDF file and convert to its internal code units (see normalization above).  Below is an example of a Python script to generate a solar wind file.

The repository includes a Jupyter Notebook, ``analysis/notebooks/Tutorial/Kaiju-IdealSW.ipynb`` that provides an interactive way to create your own solar wind file.  You can also modify the python code below to create a solar wind file.

.. code-block::

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

Visualization
-------------

Data can be read into VisIt or Paraview using the "kaiju/scripts/genXDMF.py" script which will create an XDMF file which either of those viewers can natively read.

Alternatively, you can try "kaiju/scripts/msphpic.py" script which uses the Python post-processing routines to generate a multi-panel plot of a magnetosphere run.  It uses Python libraries that you may or may not have, I don't know your life.

Creating gamera videos in parallel via slurm job submission and gamsphVid.py .

.. code-block::

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

Modify and launch the above script using "sbatch ABOVEFILE.s"

.. code-block::

   #sbatch array#

Specifies how many parallel jobs you want (36 here)

.. code-block::

   export TS="60"
   export TE="1480"
   export DT="160"

are the time domain in which you want to create movies.  TS is the start time (in minutes) of your run starting from T=0.  TE is time end (in minutes).  DT is the time (in seconds) between slices.  DT will attempt to find the cloest timestep, so if outputs are 60s but you input DT=15s, then you get duplicated frames.

.. code-block::

   source ~/.bashrc
   setdav

Set up your environment variables.  These are Kareem specific, so modify to match whatever your Casper environment is.

MPI Magnetosphere
=================

Running a basic Gamera case with MPI support requires three things:


#. Building the MPI version of the gamera executable (See [[Building_Gamera_with_MPI]])
#. Modifying case XML to supply additional MPI decomposition information
#. Modifying the submission script to request multiple nodes and use mpirun

Modifying Case XML
==================

Three additional lines are required in the case XML file when running with MPI decomposition. And these lines are ignored by the non-MPI version of Gamera, so you can safely leave them in the XML (if you want to) when not using MPI.

In the "\ :raw-html-m2r:`<Gamera>`\ " section of the XML, one line is required for each dimension that tells how many regions that dimension is decomposed into, and whether that dimension is periodic. Here is an example where the case is decomposed into 4 regions each along the I and J axes, and 2 regions along the K axis. The I and J axes are not periodic, but the K axis is.

.. code-block:: shell

   <Gamera>
   ...
      <iPdir N="4" bcPeriodic="F"/>
      <jPdir N="4" bcPeriodic="F"/>
      <kPdir N="2" bcPeriodic="T"/>
   ...
   </Gamera>

 (Corrections: Only k=1 works with the magnetosphere because the ringav code needs all the k direction data on one rank. Change to kPdir N="1". Accordingly the PBS job script should use #PBS -l select=8:ncpus=36:mpiprocs=2:ompthreads=18+1:ncpus=36:mpiprocs=1:ompthreads=36 with the part up to the plus sign allocating for gamera proper, following the plus sign allocates for voltron.)

Modifying Job Submission Script
===============================

For information about creating job submission scripts, check the `Gamerasphere <Gamerasphere>`_ article. Assuming that you have a job submission script suitable for running serial jobs, here is how to modify it to run MPI gamera. These examples will be designed to work with Cheyenne, but can be adapted to most clusters.

First, you need to request an appropriate number of nodes and tell the job submission system how many MPI ranks should be created per node. Continuing the example above, this case is decomposed 4 times along the I and J dimensions, and twice along the K dimension. So this case will need a total of 4*4*2=32 MPI ranks. For this example we will assign one MPI rank to each physical processor/socket, which provides a reasonable balance between performance and cost. Each of Cheyenne's compute nodes has two processors/sockets, so this means that each compute node will receive two MPI ranks. The original, serial, resource request line from the job submission script looked like this:

.. code-block:: shell

   #PBS -l select=1:ncpus=72:ompthreads=72

That line requests all 72 cpus on a single compute node, which is perfect for a single process. We want to create a total of 32 processes (one per MPI rank) spread across 16 compute nodes (2 MPI ranks per node), with each MPI rank getting half of the compute node's compute resources. That looks like this:

.. code-block:: shell

   #PBS -l select=16:ncpus=36:mpiprocs=2:ompthreads=36

The only other line we need to change is the one that calls the executable and starts the simulation. In the serial case that looked like this:

.. code-block:: shell

   ${EXE} ${RUNID}.xml ${JNUM} > ${RUNID}.${JNUM}.out

That command literally calls the executable and passes it the input XML file. Instead, we now need to use a helper application called mpirun, which will call our executable for us:

.. code-block:: shell

   mpirun ${EXE} ${RUNID}.xml ${JNUM} > ${RUNID}.${JNUM}.out

Using MPT
=========

If you are running with the MPT mpi library, the submission script will require some additional modifications, described in a dedicated wiki page `Running with MPT <runningWithMPT>`_.
