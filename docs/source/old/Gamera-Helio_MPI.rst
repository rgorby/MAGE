
Build/configure
===============

In addition to (instead of) instructions for the `Basic MPI Gamera <https://bitbucket.org/aplkaiju/kaiju/wiki/Building_Gamera_with_MPI>`_ do the following in your build directory to configure/build:

.. code-block:: shell

   cmake ..  (create configuration with default parameters)
   ccmake .. (ENABLE_MPI ON; GAMIC -> ...gamera/ICs/wsa.F90): [c], [g]
   make -j4 gamera_mpi.x

Then, in your run directory place your heliogrid.h5 and innerbc.h5 files and execute. An example xml files is as follows:

.. code-block:: shell

   <?xml version="1.0"?>
   <!-- Helio test -->
   <Gamera>
       <sim runid="wsa" doH5g="T" H5Grid="heliogrid.h5" icType="user"/>
       <time tFin="270."/>
       <output dtOut="10" tsOut="50.0" timer="F"/>
       <physics doMHD="T" gamma="1.5"/>
       <prob Tsolar = "25.38"/>
       <restart resFile = "wsa.Res.00008.h5" dtRes="1000." doRes="F"/>
       <iPdir N="4" bcPeriodic="F"/>
       <jPdir N="2" bcPeriodic="F"/>
       <kPdir N="4" bcPeriodic="T"/>
   </Gamera>

an example pbs script is as follows:

.. code-block:: shell

   #!/bin/bash
   #PBS -A UJHB0015
   #PBS -N heliompi
   #PBS -j oe
   #PBS -q regular
   #PBS -l walltime=01:00:00
   #PBS -l select=16:ncpus=36:mpiprocs=2:ompthreads=36
   #PBS -m ae
   #PBS -M slava.merkin@jhuapl.edu

   #Example usage

   export EXE="./gamera_mpi.x"
   export RUNID="wsa"

   source ~/.bashrc

   module list
   hostname
   date
   #export OMP_NUM_THREADS=36
   export KMP_STACKSIZE=128M
   export JNUM=${PBS_ARRAY_INDEX:-0}
   echo "Running $EXE"
   mpirun ${EXE} ${RUNID}.xml ${JNUM} > ${RUNID}.${JNUM}.out
   date
