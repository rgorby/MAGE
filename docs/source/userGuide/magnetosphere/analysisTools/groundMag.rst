
Ground Magnetometer Calculations
================================

Introduction
------------

Comparison of magnetosphere model results to ground magnetometer measurements is a common technique for validating simulations and analyzing the results. In the MAGE software, the ``calcdb.x`` program is used to calculate the magnetic field perturbations on a grid on the Earth's surface using the `Biot-Savart Law <https://en.wikipedia.org/wiki/Biot%E2%80%93Savart_law>`_, and the ionospheric, field-aligned, and magnetospheric current systems extracted from the MAGE simulation results.

This page provides an overview of how to set up and run these calculations. The :doc:`SuperMage Tools <superMAGE>` page provides instructions for conducting comparisons between these model results and data obtained from the `SuperMAG <https://supermag.jhuapl.edu/>`_ collection of ground magnetometer data.

A simple example
----------------

Preparing to run ``calcdb.x``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The program ``calcdb.x`` requires an XML file as input. The XML file provides
details on the MAGE simulation results to be used in the calculation of the
ground magnetic field perturbations, and the output requirements for the
calculation. This XML file is passed to the ``calcdb.x`` program as a
command-line argument.

Preparing the XML file
~~~~~~~~~~~~~~~~~~~~~~

Assume we have completed a simulation of the magnetosphere using the MPI MAGE code, and all results are in the current directory. The run ID is ``msphere``. The simulation results are in the HDF5 files ``msphere_*.gam.h5``. We now want to compute the ground magnetic field perturbations corresponding to this model output. Use the following specifications for the calculation:

* Start the computation at 0 simulated seconds after the start of the
  simulation results, end at 3600 simulated seconds (1 simulated hour) after
  the start of the simulation results, and provide ground magnetic field
  perturbation values at an interval of 60 simulated seconds (1 simulated
  minute).

* The spatial grid on the ground for the output will have 360 latitude bins
  and 720 longitude bins (0.5 x 0.5 degree/grid cell).

Here is a sample XML input file, followed by explanations of the individual
XML elements.

.. code-block:: XML

    <?xml version="1.0"?>
    <Kaiju>
        <Chimp>
            <sim runid="msphere"/>
            <time T0="0.0" dt="60.0" tFin="7200.0"/>
            <fields ebfile="msphere" grType="LFM" doJ="T" isMPI="T"/>
            <parallel Ri="{{ Ri }}" Rj="{{ Rj }}" Rk="{{ Rk }}"/>
        </Chimp>
    </Kaiju>

The XML elements and their attributes are described below. Note that *all* XML
attribute values are entered as strings, in ``"double quotes"``. Defaults are
supplied in ``calcdb.x`` for all values set in the XML file, so technically
this file is optional. In practice, you will want to create your own input XML
file for ``calcdb.x``\ , since the defaults for items like ``ebfile`` are not
usually appropriate for a user run.

* ``<Kaiju>`` (required): This outermost element is required.

* ``<Chimp>`` (required): This inner element is required.

* ``<sim>`` (optional): Specify identifying information for this computation.

    * ``runid`` (optional, default ``"Sim"``): String specifying an identifier
      for this run of ``calcdb.x``. A best practice is to use the ``runid`` in
      the name of the XML file.

* ``<time>`` (optional): Specify time range and interval for magnetic field
  calculation.

    * ``T0`` (optional, default ``"0.0"``): Start time (simulated seconds) for
      ground magnetic field calculation, relative to start of simulation
      results used as input.
    * ``dt`` (optional, default ``"1.0"``): Time interval and output cadence
      (simulated seconds) for ground magnetic field calculation.
    * ``tFin`` (optional, default ``"60.0"``): Stop time (simulated seconds)
      for ground magnetic field calculation, relative to start of simulation
      results used as input.

* ``<fields>`` (required): Describes the input data from a MAGE model run.

    * ``ebfile`` (optional, default ``"ebdata.h5"``): Path to HDF5 file
      containing the electric and magnetic fields computed by a MAGE model
      run.
    * ``grType`` (optional, default ``"EGG"``): String specifying grid type
      used by the MAGE output file. Valid values are ``"EGG"``\ , ``"LFM"``\ ,
      ``"SPH"``. If the string is not one of the supported grid types, the
      default value (\ ``"EGG"``\ ) is used, and a warning message is printed.
    * ``doJ`` (required, must be ``"T"``\ ): If ``"T"``\ , compute currents
      from the MAGE model results.

Once created, this XML file can be used to run the ground magnetic field
perturbation computation as follows:

.. code-block:: bash

    calcdb.x storm1.xml

The output should look something like this:

.. code-block:: bash

    goApe: Create all data structures
    Reading input deck from storm1.xml
    ----------------------------
    Initializing model ...
     KAIJU/CHIMP/units/uid                   :         Earth (DEFAULT) 
     KAIJU/CHIMP/domain/rClosed              :   3.50000E+00 (DEFAULT) 
    ------------
    CHIMP Units
    inTScl =    47.0034113607501     
    inBScl =    3.74189828029921     
    inVScl =   3.335640951981520E-006
    ------------
   ...<snip>...
   <Reading eb from msphere/Step#36>
    Writing Step#59
    Writing Step#60
   UT = 2016-08-09 10:00:00
        SMU =      263.257 [nT]
            @ Lat/Lon:       70.744      62.555 [deg]
        SML =     -219.196 [nT]
            @ Lat/Lon:       66.950     202.808 [deg]
        SME =      482.453 [nT]
        SMR =      -15.241 [nT]
            SMR-00/06/12/18 =       -9.808     -10.358     -22.392     -18.404 [nT]
          T =     3600.000 [Seconds]
      kDBps = Infinity

This command takes about 10 minutes on an i9 Mac running Monterey. When this
command completes, the file ``calcdb_storm1.out`` will contain the terminal
output from ``calcdb.x`` during the run. The computed values of the ground
magnetic field perturbation will be in the file ``storm1.deltab.h5``

Other XML elements and attributes
---------------------------------

* ``<grid>`` (optional): Options to set the grid on the ground used in
  calcdb.x, see calcdbio.F90

    * ``doH5g`` (optional, default ``"false"``\ ):  Set to ``"true"`` to use a
      grid specified in an external h5 file. If ``"false"``\ , will use a
      cartesian grid in latitude, longitude, and altitude specified by Nlat,
      Nlon,Nz
    * ``H5Grid`` (optional, default ``"grid.h5"``\ ): String specifying an the
      name of the h5 file where the grid is specified. Used if
      ``doH5g="true"``.
    * ``Nlat`` (optional, default ``"45"``\ ): Number of latitudinal cells.
    * ``Nlon`` (optional, default ``"90``\ ): Number of longitudinal cells. 
    * ``Nz`` (optional, default ``"2"``\ ): Number of cells altitude or height
      above ground
    * ``doGEO`` (optional, default ``"false"``\ ):  Set to ``"true"`` to use
      geographic coordinate system on the ground. If set to ``"false"``\ ,
      will use the SM coordinate system used in magnetosphere runs.
    * ``dzGG`` (optional, default ``"10.0"``\ ):  Height spacing [in km] of
      grid.
    * ``z0`` (optional, default ``"0.0"``\ ):  Starting height above ground
      [in km] for grid calculation.

* ``<calcdb>`` (optional): Optional setting for calcdb.x, see calcdbutils.F90

    * ``rMax`` (optional, default ``"30"``\ ): Number of latitudinal cells.
    * ``doCorot`` (optional, default ``"false"``\ ):  Set to ``"true"`` to use
      the corotation potential in the calculation.
    * ``doHall`` (optional, default ``"true"``\ ):  Set to ``"true"`` to
      include Hall currents in calculation.
    * ``doPed`` (optional, default ``"true"``\ ):  Set to ``"true"`` to
      include Pedersen currents in calculation.

* ``<parintime>`` (optional): Options to run multiple jobs of calcdb.x to increase calculation speed


  * ``NumB`` (optional, default ``"0"``\ ): Number of segments into which the data will be split for parallel computation. Must equal the number of threads requested in the PBS job script.

* 
  ``<fields>`` (required): Describes the input data from a MAGE model run.


  * ``doEBFix`` (optional, default ``"false"``\ ): Set to ``"true"`` to "clean" the electric field E so that the dot product of the electric and magnetic fields is 0. See ``ebinterp.F90``.
  * ``doMHD`` (optional, default ``"false"``\ ): Set to ``"true"`` to pass the full set of magnetohydrodynamic variables to CHIMP, rather than just electric and magnetic fields. See ``ebtypes.F90``.
  * ``isMPI`` (optional, default ``"false"``\ ): Set to ``"true"`` is the MAGE results file was generated with an MPI version of the model. See ``eblCstd.F90``.

* 
  ``<interp>`` (optional): Options related to interpolation


  * ``wgt`` (optional, default ``"TSC"``\ ): Sets 1D interpolation type. Valid values are ``"TSC"`` (1D triangular shaped cloud), ``"LIN"`` (linear),\ ``"QUAD"`` (parabolic). See ``starter.F90``.

* 
  ``<output>`` (optional): Options related to driver output


  * ``dtOut`` (optional, default ``"10.0"``\ ): Output cadence
  * ``timer`` (optional, default ``"false"``\ ): Set to ``"true"`` to turn time flags on See ``starter.F90``\ , line 139.
  * ``tsOut`` (optional, default ``"10"``\ ): Cadence to output diagnostics to run-log file See ``starter.F90``.
  * ``doFat`` (optional, default ``"false"``\ ): Set to ``"true"`` to include spherical vector components of magnetic field perturbations and currents.

* 
  ``<parallel>`` (optional): Options if ebfile was generated using an MPI version of the code (read if fields/doMPI is set to ``"true"``\ , file name in form of ebfile_Ri_Rj_Rk_i_j_k.gam.h5)


  * ``Ri`` (optional, default ``"1"``\ ): Number of ranks used in decomposition  of ``"i"`` dimension See iotable.F90.
  * ``Rj`` (optional, default ``"1"``\ ): Number of ranks used in decomposition  of ``"j"`` dimension See iotable.F90.
  * ``Rk`` (optional, default ``"1"``\ ): Number of ranks used in decomposition  of ``"k"`` dimension See iotable.F90.
  * ``doOldNaming`` (optional, default ``"false"``\ ): Allow for backward compatibility for MHD files generated with the now deprecated naming convention See ``chmpdefs.F90``.

* 
  ``<units>`` (optional): Name of units system used in the model run.


  * ``uID`` (optional, default ``"Earth"``\ ): See chmpunits.F90 line 148. Valid values are ``"EARTH"``\ , ``"EARTHCODE"``\ , ``"JUPITER"``\ , ``"JUPITERCODE"``\ , ``"SATURN"``\ , ``"SATURNCODE"``\ , ``"HELIO"``\ , ``"LFM"``\ , ``"LFMJUPITER"``.

Making the calculation go faster
--------------------------------

The computation can take a long time, so submitting this computation as a PBS job on an HPC system is usually faster, to take advantage of threads. In this case, the XML input file needs a few extra items:

.. code-block:: XML

   <?xml version="1.0"?>
   <Kaiju>
      <Chimp>
         <sim runid="storm1"/>
         <time T0="0.0" dt="60.0" tFin="39600.0"/>
         <fields doMHD="T" ebfile="msphere" grType="LFM" isMPI="T" doJ="T"/>
         <parallel Ri="4" Rj="4" Rk="1"/>
         <units uid="EARTH"/>
         <interp wgt="TSC"/>
         <grid NLat="360" NLon="720" Nz="1" doGEO="T" dzGG="10.0" z0="0.0"/>
         <output timer="F" doFat="F"/>
         <parintime NumB="32"/>
      </Chimp>
   </Kaiju>

The new item is:


* ``<parintime>`` (optional): Options to run multiple jobs of calcdb.x to increase calculation speed

  * ``NumB`` (optional, default ``"0"``\ ): Number of segments into which the data will be split for parallel computation. Must equal the number of threads requested in the PBS job script.

For example, on ``cheyenne``\ , the PBS script ``storm1.pbs`` might look something like this:

.. code-block:: bash

   #!/bin/bash
   #PBS -A P28100045
   #PBS -N storm1
   #PBS -j oe
   #PBS -q regular
   #PBS -l walltime=12:00:00
   #PBS -l select=1:ncpus=72:ompthreads=72

   #Example usage
   #qsub -v KAIJUEXE="./calcdb.x" -J 1-32 -N storm1 storm1.pbs

   # Define the kaiju installation location.                                       
   # NOTE: You should set this variable to the path to your kaiju directory.       
   export KAIJU_INSTALL_DIR=$HOME/cgs/kaiju/development/kaiju

   # Set kaiju-related environment variables.                                      
   source $KAIJU_INSTALL_DIR/scripts/setupEnvironment.sh

   # Add the kaiju binary directory to the command path.                           
   export EXE=${KAIJUEXE:-"./calcdb.x"}
   export RUNID=${PBS_JOBNAME}

   # Load the required modules for MPI kaiju.
   # Replace this list with the modules used to build your kaiju installation.
   module purge
   module restore kaiju2021

   module list
   hostname
   date
   export OMP_NUM_THREADS=72
   export KMP_STACKSIZE=128M
   export JNUM=${PBS_ARRAY_INDEX:-0}
   echo "Running $EXE"
   echo "RUNID is $RUNID"
   ${EXE} ${RUNID}.xml ${JNUM} > ${RUNID}.${JNUM}.out
   date

This job can be submitted with the command

.. code-block:: bash

   qsub -v KAIJUEXE="./calcdb.x" -J 1-32 -N storm1 storm1.pbs

Note that the value of the number of jobs requested must match the value of the ``NumB`` attribute of the ``<parintime/>`` element from the XML file.

When this job completes, files of the form ``calcdb_storm1.#.out`` which will contain the terminal output from ``calcdb.x`` during the run for each job. 

Finishing up
------------

After your batch completes, you'll need to do one more step before you can
plot your results or use :doc:`SuperMage Tools <superMAGE>` to make
comparisons with ground magnetometers.

Concatenating multiple HDF5 files from MPI processing
-----------------------------------------------------

The parallel processing results in multiple output files. The ``dbCat.py`` script in the ``scripts/postproc`` directory of the ``kaiju`` repository will concatenate them into one file. The computed values of the ground magnetic field perturbation will be in the file ``storm1.deltab.h5``.

.. code-block:: bash

   dbCat.py -runid storm1
