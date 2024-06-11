.. role:: raw-html-m2r(raw)
   :format: html


MPI Differences
===============

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

 (Corrections: Only k=1 works with the magnetosphere because the ringav code needs all the k direction data on one rank. Change to kPdir N="1". Accordingly the PBS job script should use #PBS -l select=16:ncpus=36:mpiprocs=2:ompthreads=18+1:ncpus=36:mpiprocs=1:ompthreads=36 with the part up to the plus sign allocating for gamera proper, following the plus sign allocates for voltron.)

Modifying Job Submission Script
===============================

For information about creating job submission scripts, check the [[Gamerasphere]] article. Assuming that you have a job submission script suitable for running serial jobs, here is how to modify it to run MPI gamera. These examples will be designed to work with Cheyenne, but can be adapted to most clusters.

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

If you are running with the MPT mpi library, the submission script will require some additional modifications, described in a dedicated wiki page [[Running with MPT]].
