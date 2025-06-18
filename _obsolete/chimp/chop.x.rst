Extract a 3D subdomain from GAMERA output (chop.x)
==================================================

Introduction
------------

chop.x extracts a 3D portion of the domain from GAMERA output on either the
MAGE grid or interpolated onto a cartesian or a spherical grid and perform
additional calculations such as field line tracing. Output is to an hdf5 file
that can be visualized and analyzed similarly to GAMERA output files.

Example XML file
----------------

.. code-block::

   <?xml version="1.0" ?>
   <Kaiju>
     <Chimp>
      <sim runid="chopRTP"/>
      <time T0="26100.0" dt="30.0" tFin="26520.1"/>
      <fields doMHD="T" ebfile="msphere" grType="LFM" isMPI="T"/>
      <parallel Ri="8" Rj="8" Rk="1"/>
      <domain dtype="LFM" xSun="30.0" xTail="-50.0" yzMax="30.0"/>
      <units uid="EARTH"/>
      <chop grType="RTP" x1Max="20" x1Min="2" x2Max="180" x2Min="0" x3Max="360" x3Min="0" Nx1="225" Nx2="120" Nx3="240"/>
      <interp wgt="TSC"/>
      <tracer epsds="0.20"/>
      <output doTrc="T"/>
     </Chimp>
   </Kaiju>

Parameter Descriptions
----------------------


``<sim>`` (optional): Specify identifying information for this computation.

``runid`` (optional, default ``"Sim"``): String specifying an identifier
for this run of ``chop.x``. A best practice is to use the ``runid`` in the
name of the XML file.

``<time>`` (optional): Specify time range and interval for magnetic field
calculation.

``T0`` (optional, default ``"0.0"``): Start time (simulated seconds) for
ground magnetic field calculation, relative to start of simulation results
used as input.

``dt`` (optional, default ``"1.0"``): Time interval and output cadence
(simulated seconds) for ground magnetic field calculation.

``tFin`` (optional, default ``"60.0"``): Stop time (simulated seconds) for
ground magnetic field calculation, relative to start of simulation results
used as input.

``<fields>`` (required): Describes the input data from a MAGE model run.

``ebfile`` (optional, default ``"ebdata.h5"``): Path to HDF5 file
containing the electric and magnetic fields computed by a MAGE model run.

``grType`` (optional, default ``"EGG"``): String specifying grid type used
by the MAGE output file. Valid values are ``"EGG"``, ``"LFM"``, ``"SPH"``.
If the string is not one of the supported grid types, the default value
(``"EGG"``) is used, and a warning message is printed.

``doEBFix`` (optional, default ``"false"``): Set to ``"true"`` to "clean"
the electric field E so that the dot product of the electric and magnetic
fields is 0. See ``ebinterp.F90``.

``doMHD`` (optional, default ``"false"``): Set to ``"true"`` to pass the
full set of magnetohydrodynamic variables to CHIMP, rather than just electric
and magnetic fields. Includes velocity vector, density and pressure in the
output file. See ``ebtypes.F90``.

``isMPI`` (optional, default ``"false"``): Set to ``"true"`` is the MAGE
results file was generated with an MPI version of the model. See
``eblCstd.F90``.

``<domain>`` (optional): Options for the problem domain

``dtype`` (optional, default ``"SPH"``): Domain over which to perform
CHIMP calculations, separate from grid, enables the user to perform
calculation on a subset of the grid to reduce computation where it is not
needed - See ``gridloc.F90``. Valid values are ``"SPH"``, ``"LFM"``,
``"LFMCYL"``, ``"MAGE"``, ``"EGG"``, ``"ELL"``.

``rClosed`` (optional, default set by choice of ``units/uid``): Radial value
for field line endpoint to reach to be considered closed - See
``chmpunits.F90``.

``rmax`` (optional, default computed): Maximum radius of Domain region -
See ``gridloc.F90``.

``rmin`` (optional, default computed): Minimum radius of Domain region -
See ``gridloc.F90``.

``xSun`` (optional,default 20.0): if dType is "LFM" or "MAGE", the Domain
region includes all i-shells which have distances along the Earth-Sun line is
less than this value (in Re)

``xTail`` (optional,default -100.0): if dType is "LFM" or "MAGE", the Domain
region includes cells in the magnetotail up until this value (in Re)

``yzMax`` (optional,default 40.0): if dType is "LFM" or "MAGE", the Domain
region includes cells with Y and Z coordinates between +/- yzMax (in Re)

``<output>`` (optional): Options related to driver output

``timer`` (optional, default ``"false"``): Set to ``"true"`` to turn time
flags on See ``starter.F90``.

``tsOut`` (optional, default ``"10"``): Cadence to output diagnostics to
run-log file See ``starter.F90``.

``doEQProj`` (optional, default ``"false"``): Set to ``.true.`` to include
equatorial variables, projected down to magnetic equator along field line from
cell location (i.e. Xeq,Yeq, if field line is open or closed etc) See
``chmpdefs.F90``.

``doSlim`` (optional, default ``"false"``):  Set to ``"true"`` to remove
vector electric field and current data from slice.x output See
``chmpdefs.F90``.

``doTrc`` (optional, default ``"false"``): Similar to doEQProj, used in
chop.x See ``chmpdefs.F90``.

``<parallel>`` (optional): Options if ebfile was generated using an MPI
version of the code (read if fields/doMPI is set to ``"true"``, file name in
form of ebfile_Ri_Rj_Rk_i_j_k.gam.h5)

``Ri`` (optional, default ``"1"``): Number of ranks used in decomposition of
``"i"`` dimension See iotable.F90.

``Rj`` (optional, default ``"1"``): Number of ranks used in decomposition of
``"j"`` dimension See iotable.F90.

``Rk`` (optional, default ``"1"``): Number of ranks used in decomposition of
``"k"`` dimension See iotable.F90.

``doOldNaming`` (optional, default ``"false"``): Allow for backward
compatibility for MHD files generated with the now deprecated naming
convention See ``chmpdefs.F90``.

``<tracer>`` (optional): Options related to field line tracing performed by
CHIMP

``epsds`` (optional, default ``"1.0e-2"``): Tolerance for field line tracing
computations See chmpdefs.F90.

``<units>`` (optional): Name of units system used in the model run.

``uID`` (optional, default ``"Earth"``): See chmpunits.F90. Valid values
are ``"EARTH"``, ``"EARTHCODE"``, ``"JUPITER"``, ``"JUPITERCODE"``,
``"SATURN"``, ``"SATURNCODE"``, ``"HELIO"``, ``"LFM"``, ``"LFMJUPITER"``.

``<interp>`` (optional): Options related to interpolation

``wgt`` (optional, default ``"TSC"``): Sets 1D interpolation type. Valid
values are ``"TSC"`` (1D triangular shaped cloud), ``"LIN"`` (linear),
``"QUAD"`` (parabolic). See ``starter.F90``.

``<chop>`` (optional): Options specific to chop.x driver, see chopio.F90

``grType`` (optional, default ``"XYZ"``): String specifying an identifier
for the grid to do 3D data extraction. Valid Values are ``"XYZ"`` (cartesian),
``"RTP"`` (spherical), ``"LFM"`` (MAGE grid).

``Nx1`` (optional, default ``"64"``): Number of cells in X or R depending on
grid specified. Not used if grType is set to ``"LFM"``.

``Nx2`` (optional, default ``"64"``): Number of cells in Y or Theta
depending on grid specified. Not used if grType is set to ``"LFM"``.

``Nx3`` (optional, default ``"64"``): Number of cells in Z or Phi depending
on grid specified. Not used if grType is set to ``"LFM"``.

``x1Max`` (optional, default ``"10.0"``): Maximum value of the X1 dimension
used to initialize the chop grid. If ``grTyp="LFM"``, x1Max is used similar to
domain/xSun

``x1Min`` (optional, default ``"-x1Max"``): Minimum value of the X1
dimension used to initialize the chop grid. Not used if grType is set to
``"LFM"``.

``x2Max`` (optional, default ``"10.0"``): Maximum value of the X2 dimension
used to initialize the chop grid. Not used if grType is set to ``"LFM"``.

``x2Min`` (optional, default ``"-x2Max"``): Minimum value of the X2
dimension used to initialize the chop grid. Not used if grType is set to
``"LFM"``.

``x3Max`` (optional, default ``"10.0"``): Maximum value of the X3 dimension
used to initialize the chop grid. Not used if grType is set to ``"LFM"``.

``x3Min`` (optional, default ``"-x3Max"``): Minimum value of the X3
dimension used to initialize the chop grid. Not used if grType is set to
``"LFM"``.

Run Script example
------------------

An example pbs script on ``cheyenne``, RunChop.pbs to submit a chop job:

.. code-block:: bash

   #!/bin/bash
   #PBS -A P28100045
   #PBS -N chopRTP
   #PBS -j oe
   #PBS -q regular
   #PBS -l walltime=12:00:00
   #PBS -l select=1:ncpus=72:ompthreads=72
   export EXE=${chop}
   export RUNID=${PBS_JOBNAME}

   #Replace this with your module set
   module purge
   module restore kaiju

   module list
   hostname
   date
   export OMP_NUM_THREADS=72
   export KMP_STACKSIZE=128M
   echo "Running $EXE"
   ./chop.x ${RUNID}.xml > ${RUNID}.out
   date

This job can be submitted with the command

.. code-block:: bash

   qsub RunChop.pbs
