CHIMP Quick Start Guide
=======================

Introduction
------------

This page provides instructions to help the user get up and running quickly
with `CHIMP <https://cgs.jhuapl.edu/Models/chimp.php>`_.

The tools in CHIMP are:

    * ``push.x`` - Computes the trajectory of test particles within the domain
      of the MAGE model run.
    * ``slice.x`` - Extracts 4-D slices from the MAGE results.
    * ``trace.x`` - Traces field lines in the MAGE results.

These instructions assume you are running on an HPC system that uses PBS, and
that you have already made a successful ``kaiju`` model run in your current
working directory. Each tool requires an XML file to run.

Running ``push.x``
------------------

Create the XML file
^^^^^^^^^^^^^^^^^^^

Create an XML file describing the input parameters for ``push.x``. It should
look something like this on ``derecho``:

.. code-block:: xml

   <?xml version="1.0" ?>
   <Kaiju>
      <Chimp>
         <sim runid="geospace"/>
         <time T0="36000.0" dt="0.5" tFin="43200.0"/>
         <fields doMHD="F" ebfile="geospace" grType="LFM" isMPI="T"/>
         <parallel Ri="4" Rj="4" Rk="1"/>
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

Notes: Here ``T0`` is start time for the ``push.x`` run taken as seconds from
LFM's T=0s. ``tFin`` is the stop time for ``push.x``. For other parameters
then see the :doc:`Chimp XML </userGuide/chimp/chimpXML>` page.

Create the PBS file
^^^^^^^^^^^^^^^^^^^

Your PBS script to run ``push.x`` should look something like this:

.. code-block:: bash

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

Run the tool
^^^^^^^^^^^^

Create your run folder and copy chimp's executable in.

.. code-block:: bash

   cd $SCRATCH
   mkdir TESTRUN
   cd TESTRUN
   cp $KAIJUDIR/build/bin/push.x ./

An example command to submit two parallel jobs:

.. code-block:: bash

   qsub -v CHIMPEXE="./push.x" -J 1-2 -N Oxyflow RunCHIMP.pbs


Notes:

a. Need to add dot slash (``./``) in front of the executable name so
that it can be found, or you can add ``./``, namely the current working
directory to your path.

b. ``-J 1-2`` means the same job is submitted in parallel two times.

c. ``-N Oxyflow`` specifies the run name is ``Oxyflow``, which has to be
consistent with the xml file name ``Oxyflow.xml`` (don't include ``.xml`` when
running ``qsub``). This name will also overwrite whatever in the PBS file.
Example XML and PBS files can be found in
:doc:`here </userGuide/chimp/chimpXML>`.

The job number is used to seed the generation of test particles within each
run. This allows the code to be reproducible. If you want to run additional
test particles with the same xml file, you need to make sure the job numbers
are not the same i.e. if your first set of runs was submitted with the command
above, you can run two additional, unique jobs with the command:

.. code-block:: bash

   qsub -v CHIMPEXE="./push.x" -J 3-4 -N Oxyflow RunCHIMP.pbs

6. Now we run chimp. The magic happens so fast that you'll see the run
complete by the time you type in ls.

.. code-block:: bash

   qsub -v CHIMPEXE="./push.x" -J 1-2 -N TESTRUNNAME RunCHIMP.pbs


Note that the pbs file makes chimp look for ``TESTRUNNAME.xml`` so you'll want
to make sure those are consistent.

h5part Variable Descriptions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The output file from a CHIMP run will be look something like
``runid.000000.h5part`` and contain the variables:

``x``, ``y``, ``z``: particle position coordinates

``K``: particle energy [keV]

``isIn``: 0,1 depending on whether the particle is alive or not (left the
simulation domain) and whether it's been born yet or not

``id``: particle ID #

``xeq``, ``yeq``, ``Teq``: position where the particle last crossed the Z=0
plane, and at what time or instantaneous projection depending on doTrc in xml.

``Keq``/``ebKeq``: energy at the equator in lab and ExB frame

``Aeq``: equatorial pitch angle

``OCb``: is ``0``/``1``/``2`` based on the topology of the field line at the
particle position: IMF/open/closed. this will be all zeroes if you didn't tell
it to do field-line tracing (expensive) in the XML file

``T0p`` : is the birthday of the particle, when it was/will be born

Running ``slice.x``
-------------------

Running ``trace.x``
-------------------
