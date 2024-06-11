
CHIMP
=====

1. Compile chimp after gamera has been compiled successfully.

.. code-block::

   #!shell
   module purge
   module restore kaiju
   cd $KAIJUDIR/build
   cmake ..
   make gamera
   make chimp


The executables can be found in $KAIJUDIR/build/bin

.. code-block::

   #!shell
   ls $KAIJUDIR/build/bin
   gamera.x*  project.x*  psd.x*  push.x*  slice.x*


Notes: a) "module restore kaiju" is to load the 8 critical modules
before compiling (see
https://bitbucket.org/aplkaiju/kaiju/wiki/Quick_Start). b) Cheyenne
users need to make sure ncar_pylib is activated, and may sometimes need
to deactivate and then re-activate.

2. Run chimp.

First check the prerequisites are all present in the work directory,
including: a) msphere files from GAMERA output; b) xml file specifying
the chimp parameters; c) pbs file for submitting jobs if using Cheyenne
etc.; d) executable "push.x". An example command to submit two parallel
jobs:

.. code-block::

   #!shell
   qsub -v CHIMPEXE="./push.x" -J 1-2 -N Oxyflow RunCHIMP.pbs


Notes: a) Need to add dot slash (./) in front of the executable name so
that it can be found, or you can add "./", namely the current work
directory to your path. b) "-J 1-2" means the same job is submitted in
parallel two times. c) "-N Oxyflow" specifies the run name is "Oxyflow",
which has to be consistent with the xml file name Oxyflow.xml (don't
include .xml when qsub). This name will also overwrite whatever in the
pbs file. Example xml and pbs files can be found in
"https://bitbucket.org/aplkaiju/kaiju/wiki/ChimpXML".

The job number is used to seed the generation of test particles within
each run. This allows the code to be reproducible. If you want to run
additional test particles with the same xml file, you need to make sure
the job numbers are not the same i.e. if your first set of runs was
submitted with the command above, you can run two additional, unique
jobs with the command:

.. code-block::

   #!shell
   qsub -v CHIMPEXE="./push.x" -J 3-4 -N Oxyflow RunCHIMP.pbs



#. CHIMP-LFM # Anything written in all caps should be changed to match
   your run.


.. raw:: html

   <!-- end list -->




#. Compile CHIMP using the above instructions.
#. Convert an existing LFM run into a chimp compatible run using the
   built in conversion script. The --mhd option supposedly carries the
   rest of the mhd quantities too, but not necessarily needed. By
   default, this creates a file called ebLFM.h5. The output file name
   can be changed via ' -o OUTPUTFILENAME ' option.


.. raw:: html

   <!-- end list -->



.. code-block::

   #!shell
   cd $YOURLFMRUN
   lfm2chimp.py --mhd *mhd*hdf


Notes: LFM usually has an extra spinup file that is 50 minutes prior to
T=0s. You'll want to remove this prior to running the above command
otherwise it'll screw up what is suppose to be T=0s. Alternatively,
adjust all times to account for the extra 50 minutes.

3. Create your run folder and copy chimp's executable in.

.. code-block::

   #!shell
   cd $SCRATCH
   mkdir TESTRUN
   cd TESTRUN
   cp $KAIJUDIR/build/bin/push.x ./


4. Create an XML file. You can create one based on the information
above or modify mine. The XML file should be named TESTRUNNAME.xml .

.. code-block::

   #!xml
   <?xml version="1.0" ?>
   <Kaiju>
       <Chimp>
           <sim runid="Oxyflow"/>
           <time T0="36000.0" dt="0.5" tFin="43200.0"/>
           <fields doMHD="F" ebfile="ebLFM" grType="LFM" isMPI="F"/>
           <parallel Ri="6" Rj="12" Rk="1"/>
           <pusher epsht="0.05" imeth="FO"/>
           <tps Np="2500" species="Op"/>
           <units uid="LFM"/>
           <output doEQProj="T" dtOut="5.0" tsOut="100"/>
           <radius max="2.25" min="2.05"/>
           <phi max="360.0" min="0.0"/>
           <alpha max="180.0" min="90.0"/>
           <height max="90.0" min="15.0" doOutflow="T"/>
           <energy doEBInit="F" doLog="F" max="10.0" min="0.001"/>
           <domain rmax="25.0" rmin="2.005"/>
           <tracer epsds="0.05"/>
           <stream doStream="F" max="39600.0" min="36000.0"/>
       </Chimp>
   </Kaiju>


Notes: Here T0 is start time for the chimp run taken as seconds from
LFM's T=0s. tFin is chimp's stop time. For other parameters then see the
chimpXML wiki. If you made changes to the eb file name then you'll want
to change ebfile accordingly.

5. Now you'll want to create a submission environment for your job. For
this example, I'll create one based on Cheyenne (RunCHIMP.pbs) so you'll
want to create a similar one for your environment.

.. code-block::

   #!shell

   #!/bin/bash
   #PBS -A P28100045
   #PBS -N CHIMP
   #PBS -j oe
   #PBS -q regular
   #PBS -l walltime=12:00:00
   #PBS -l select=1:ncpus=72:ompthreads=72
   #PBS -m ae

   #Example usage
   #qsub -v CHIMPEXE="pusher.x" -J 1-5 -N RC_ep RunCHIMP.pbs

   export EXE=${CHIMPEXE:-"push.x"}
   export RUNID=${PBS_JOBNAME}

   module list
   hostname
   date
   export OMP_NUM_THREADS=72
   export KMP_STACKSIZE=128M
   export JNUM=${PBS_ARRAY_INDEX:-0}
   echo "Running $EXE"
   ${EXE} ${RUNID}.xml ${JNUM} > ${RUNID}.${JNUM}.out
   date


6. Now we run chimp. The magic happens so fast that you'll see the run
complete by the time you type in ls.

.. code-block::

   #!shell
   qsub -v CHIMPEXE="./push.x" -J 1-2 -N TESTRUNNAME RunCHIMP.pbs


Note that the pbs file makes chimp look for TESTRUNNAME.xml so you'll
want to make sure those are consistent.

h5part Variable Descriptions
----------------------------

The output file from a CHIMP run will be look something like \<run
id&gt;.000000.h5part and contain the variables:

x,y,z : particle position coordinates

K : particle energy [keV]

isIn : 0,1 depending on whether the particle is alive or not (left the
simulation domain) and whether it's been born yet or not

id : particle ID #

xeq,yeq,Teq : position where the particle last crossed the Z=0 plane,
and at what time or instantaneous projection depending on doTrc in xml

Keq/ebKeq : energy at the equator in lab and ExB frame

Aeq : equatorial pitch angle

OCb : is 0/1/2 based on the topology of the field line at the particle
position: IMF/open/closed. this will be all zeroes if you didn't tell it
to do field-line tracing (expensive) in the XML file

T0p : is the birthday of the particle, when it was/will be born
