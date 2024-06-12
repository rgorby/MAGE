Testing the kaiju software
==========================

Running a test case (serial)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Before you begin
~~~~~~~~~~~~~~~~

Before proceeding, initialize your ``kaiju`` environment by loading the required modules. For the serial version of ``kaiju`` on Pleiades, these commands should work:

.. code-block::

   ##!shell
   module purge
   module load pkgsrc/2021Q2
   module load comp-intel/2020.4.304
   module load hdf5/1.8.18_serial


For the MPI version of ``kaiju`` on Pleiades, use:

.. code-block::

   #!shell
   module purge
   module load pkgsrc/2021Q2
   module load comp-intel/2020.4.304
   module load mpi-hpe/mpt.2.23
   module load hdf5/1.8.18_mpt


Next, you need to run the setup script for the CDF software used by the ``kaiju`` software:

.. code-block::

   #!shell
   source CDF_INSTALL_DIR/bin/definitions.B


where ``CDF_INSTALL_DIR`` is the installation directory for your CDF software. The ``definitions.B`` script is for the ``bash`` shell. Setup scripts for other shells are in the same directory.

At this point, you should make sure your python environment is properly configured. For these examples, we assume the use of a ``conda``\ -based virtual environment for Python 3.8 called ``kaiju-python-3.8``\ , which should contain all of the required python modules described in LINK_TO_PYTHON_SECTION. Assuming this environment was already created, activate it with the command:

.. code-block::

   #!shell
   conda activate kaiju-python-3.8


Next, run the ``kaiju`` setup script appropriate for your shell. For example, using the ``bash`` shell on Pleiades, run:

.. code-block::

   #!shell
   source KAIJU_INSTALL_DIR/scripts/setupEnvironment.sh


where ``KAIJU_INSTALL_DIR`` is the path to the directory created when you cloned the ``kaiju`` repository. This script sets the ``KAIJUHOME`` environment variable (it will be the same as ``KAIJU_INSTALL_DIR``\ ), and adds ``kaiju``\ -specific entries in the ``PATH`` and ``PYTHONPATH`` environment variables.

Finally, some of the ``kaiju`` scripts require the ``geopack-2008`` compiled Python module. Add it to your Python module search path using the command:

.. code-block::

   #!shell
   export PYTHONPATH=$KAIJUHOME/external/geopack-2008/build/lib.linux-x86_64-3.8:$PYTHONPATH


where the ``lib.XXX`` string is for Pleiades; it may differ for your machine.

NOTE: The scripts used below are available as part of the ``kaiju`` code distribution. They were designed to illustrate the steps needed to run a simple model and examine the results. The scripts illustrate how to set up a run, and how to extract data from the result files. Feel free to use these scripts as the starting point for your own scripts using the ``kaiju`` code.

NOTE: All examples below assume use of the ``bash`` shell. Modify as needed for your preferred shell. The primary differences will 
be replacing ``export X=Y`` with ``setenv X Y`` when setting environment variables in the ``csh``\ /\ ``tcsh`` shells.

2-D field loop convection
~~~~~~~~~~~~~~~~~~~~~~~~~

The first example we will run is the ``loop2d`` model. This is a simple 2-D model illustrating the convective transport of a loop-shaped magnetic field created by a linear current. This test case is discussed at `Field Loop Test <https://www.astro.princeton.edu/~jstone/Athena/tests/field-loop/Field-loop.html>`_. This test case uses the *serial* version of the ``kaiju`` code.

Begin by adding the directory containing the serial ``kaiju`` binaries to your ``PATH`` environment variable.

.. code-block::

   #!shell
   export PATH=$KAIJUHOME/build/bin:$PATH


Next, add the directory containing the the ``loop2d`` example scripts to your ``PATH`` environment variable: 

.. code-block::

   #!shell
   export PATH=$KAIJUHOME/quickstart/loop2d:$PATH


Now create a new directory to run the ``loop2d`` example (it can be anywhere).

.. code-block::

   #!shell
   cd $HOME
   mkdir -p kaiju_test/loop2d
   cd kaiju_test/loop2d


The next step is to generate the configuration files for the ``loop2d`` model. This is done using the ``prepare_loop2d.py`` utility:

.. code-block::

   #!shell
   prepare_loop2d.py -v


This command will create 3 files in your current directory:


* ``loop2d.ini`` - An .ini-format initialization file for running ``gamera.x`` on the ``loop2d`` example.
* ``loop2d.xml`` - An XML-format initialization file for running ``gamera.x`` on the ``loop2d`` example, created from ``loop2d.ini``.
* ``loop2d.pbs`` - An PBS job script for running the ``gamera.x`` binary on the ``loop2d`` example.

NOTE: The conversion from ``.ini`` to ``.xml`` is still under development, so both initialization files are created from templates.

NOTE: The PBS script is designed for use on the Pleiades system at NASA Ames. If you are working in a non-HPC environment, the commands listed in the ``loop2d.pbs`` file can be executed manually on your command line.

Next, submit the PBS job script for execution. On Pleiades, the job is submitted using the ``qsub`` command:

.. code-block::

   #!shell
   qsub loop2d.pbs


Once the job has been accepted in the queue, the run should take about 20-30 seconds (on Pleiades).

When complete, you should see in your working directory the input files created by ``prepare_loop2d.py``\ , and the output files created by ``gamera.x``\ :

.. code-block::

   #!shell
   bash-4.2$ ls -1
   gamera.x.loop2d.out
   loop2d.gam.Res.00000.h5
   loop2d.gam.Res.XXXXX.h5
   loop2d.gam.h5
   loop2d.ini
   loop2d.o12913999
   loop2d.pbs
   loop2d.xml


NOTE: The ``gamera.x.loop2d.out`` contains the terminal output generated by the ``gamera.x`` executable. The file ``loop2d.o12913999`` (the digits will be different for your run) contains the PBS log for the job. The ``.h5`` files are the HDF5 output files created by ``gamera.x``.

We can now analyze the results of the model run. A simple example analysis can be run using the utility ``run_loop2d_checks.py``\ :

.. code-block::

   #!shell
   $ run_loop2d_checks.py -v
   Computing volume-integrated magnetic pressure.
   Volume-integrated magnetic pressure (SUM(Pb*dV), code units):
   At start: 1.3978992818223555e-07
   At end: 1.3276572728297693e-07


Your values for the volume-integrated magnetic pressure should be very close to these values. If they are significantly different, please double-check your build and run procedure. If you are unable to identify the cause of the discrepancy, please contact the ``kaiju`` team.

Finally, generate a quick-look plot illustrating the model results. For this case, the quick-look plot shows the magnetic pressure in the first and last simulation frames:

.. code-block::

   #!shell
   $ create_loop2d_quicklook.py -v


This script will create the file ``loop2d_quicklook.png``. It should look like this:


.. image:: https://bitbucket.org/repo/kMoBzBp/images/720765616-loop2d_quicklook.png
   :target: https://bitbucket.org/repo/kMoBzBp/images/720765616-loop2d_quicklook.png
   :alt: loop2d_quicklook.png


**Note for Cheyenne users** Sometimes you may get an error message of
*"No module h5py"*. If ncar_pylib has been loaded, try to
##deactivate## and then re-activate.

Geospace example (serial)
~~~~~~~~~~~~~~~~~~~~~~~~~

An additional serial test case is available, in the ``geo_serial`` model. This model examines the terrestrial magnetosphere during the period 2016-08-09T09:00:00 to 2016-08-09T11:00:00.

Begin by adding the directory containing the serial ``kaiju`` binaries to your ``PATH`` environment variable.

.. code-block::

   #!shell
   export PATH=$KAIJUHOME/build/bin:$PATH


To use the ``geo_serial`` test case:

.. code-block::

   #!shell
   export PATH=$KAIJUHOME/quickstart/geo_serial:$PATH


This is a more complex example of using the serial ``kaiju`` code, as it includes several preprocessing steps, and uses the ``voltron.x`` binary rather that the ``gamera.x`` binary.

The next step is to generate the configuration and input files for the ``geo_serial`` model. This is done using the ``prepare_geo_serial.py`` utility:

.. code-block::

   #!shell
   prepare_geo_serial.py -v


This command will create 7 files in your current directory:


* ``geo_serial.ini`` - An .ini-format initialization file for running ``voltron.x`` on the ``geo_serial`` example.
* ``geo_serial.pbs`` - An PBS job script for running ``voltron.x`` on the ``geo_serial`` example.
* ``geo_serial.xml`` - An XML-format initialization file for running ``voltron.x`` on the ``geo_serial`` example, created from ``geo_serial.ini``.
* bcwind.h5 - A HDF5 file containing solar wind data from `OMNIWeb <https://omniweb.gsfc.nasa.gov/ow.html>`_ to use as a boundary condition for the simulation.
* lfmD.h5 - A HDF5 file containing the grid to use for the simulation.
* OMNI_HRO_1MIN.txt.png - A quick-look plot of the solar wind data from OMNIWeb.
* rcmconfig.h5 - A HDF5 file containing configuration parameters for the Rice Convection Model (RCM) used by ``voltron.x``.

NOTE: The detailed steps for manually generating these files are described `here <geoQuickStart>`_.

NOTE: The conversion from ``.ini`` to ``.xml`` is still under development, so both initialization files are created from templates.

The PBS script is designed for use on the Pleiades system at NASA Ames. If you are working in a non-HPC environment, the commands listed in this file can be executed manually on your command line. On Pleiades, the job is submitted using the ``qsub`` command:

.. code-block::

   #!shell
   qsub geo_serial.pbs


Once the job has been accepted in the queue, the run should take about an hour (on Pleiades).

When complete, you should see in your working directory the output files created by ``voltron.x``\ , along with the files created by ``prepare_geo_serial.py``\ :

.. code-block::

   #!shell
   bash-4.2$ ls -1
   OMNI_HRO_1MIN.txt.png
   bcwind.h5
   geo_serial.ini
   geo_serial.o12920410
   geo_serial.pbs
   geo_serial.xml
   lfmD.h5
   msphere.RCM.Res.00000.h5
   msphere.RCM.Res.00001.h5
   msphere.RCM.Res.00002.h5
   msphere.RCM.Res.XXXXX.h5
   msphere.gam.Res.00000.h5
   msphere.gam.Res.00001.h5
   msphere.gam.Res.00002.h5
   msphere.gam.Res.XXXXX.h5
   msphere.gam.h5
   msphere.mhd2imag.Res.00000.h5
   msphere.mhd2imag.Res.00001.h5
   msphere.mhd2imag.Res.00002.h5
   msphere.mhd2imag.Res.XXXXX.h5
   msphere.mhdrcm.h5
   msphere.mix.Res.00000.h5
   msphere.mix.Res.00001.h5
   msphere.mix.Res.00002.h5
   msphere.mix.Res.XXXXX.h5
   msphere.mix.h5
   msphere.rcm.h5
   msphere.volt.Res.00000.h5
   msphere.volt.Res.00001.h5
   msphere.volt.Res.00002.h5
   msphere.volt.Res.XXXXX.h5
   msphere.volt.h5
   rcmconfig.h5
   voltron.x.geo_serial.out


The ``voltron.x.geo_serial.out`` contains the terminal output generated by the ``voltron.x`` executable. The file ``geo_serial.o12920410`` (the digits will be different) contains the PBS log for the job.

We can now analyze the results of the model run. A simple analysis can be run using the utility ``run_geo_serial_checks.py``\ :

.. code-block::

   #!shell
   $ run_geo_serial_checks.py -v --runid=msphere
   Computing volume-integrated magnetic pressure.
   Volume-integrated magnetic pressure (SUM(Pb*dV), code units):
   At start: 450951034.31847197
   At end: 742967690.5682685


Your values for the volume-integrated magnetic pressure should be very close to these values. This is actually a dummy statistic - it is not scientifically useful, but the code illustrates how to access and manipulate the results using the ``kaiju`` software.

Finally, generate a quick-look plot illustrating the model results. For this case, the quick-look plot shows several plots of different data generated by the model:

.. code-block::

   #!shell
   $ create_geo_serial_quicklook.py -v --runid=msphere


This script will create the file ``qkpic.png``. It should look like this:


.. image:: https://bitbucket.org/repo/kMoBzBp/images/4238837198-qkpic.png
   :target: https://bitbucket.org/repo/kMoBzBp/images/4238837198-qkpic.png
   :alt: qkpic.png


Heliosphere example (serial)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``helio_serial`` example is a heliospheric model, using the heliosphere-specific build of ``gamera`` (\ ``gamhelio.x``\ ).

Begin by adding the directory containing the serial ``kaiju`` binaries to your ``PATH`` environment variable.

.. code-block::

   #!shell
   export PATH=$KAIJUHOME/build/bin:$PATH


To use the ``helio_serial`` test case, add the model directory to your command path:

.. code-block::

   #!shell
   export PATH=$KAIJUHOME/quickstart/helio_serial:$PATH


The next step is to generate the configuration and input files for the ``helio_serial`` model. This is done using the ``prepare_helio_serial.py`` utility:

.. code-block::

   #!shell
   prepare_helio_serial.py -v


This command will create several files in your current directory:


* ``helio_serial.ini`` - An .ini-format initialization file for running ``gamhelio.x`` on the ``helio_serial`` example.
* ``helio_serial.pbs`` - An PBS job script for running ``gamhelio.x on the``\ helio_serial` example.
* ``helio_serial.xml`` - An XML-format initialization file for running ``gamhelio.x`` on the ``helio_serial`` example, created from ``helio_serial.ini``.
* heliogrid.h5 - A HDF5 file containing the grid to use for the simulation.
* innerbc.h5 - A HDF5 file containing the inner boundary conditions derived from the the WSA (Wang-Sheeley-Arge) model used for this example.

NOTE: The detailed steps for manually generating these files are described `here <helioQuickStart>`_.

NOTE: The conversion from ``.ini`` to ``.xml`` is still under development, so both initialization files are created from templates.

The PBS script is designed for use on the Pleiades system at NASA Ames. If you are working in a non-HPC environment, the commands listed in this file can be executed manually on your command line. On Pleiades, the job is submitted using the ``qsub`` command:

.. code-block::

   #!shell
   qsub helio_serial.pbs


Once the job has been accepted in the queue, the run should take about 30 minutes (on Pleiades).

When complete, you should see in your working directory the output files created by ``gamhelio.x``\ , along with the files created by ``prepare_helio_serial.py``\ :

.. code-block::

   #!shell
   bash-4.2$ ls -1
   gamhelio.x.helio_serial.out
   helio_serial.ini
   helio_serial.o13122765
   helio_serial.pbs
   helio_serial.xml
   heliogrid.h5
   innerbc.h5
   wsa.gam.Res.00000.h5
   wsa.gam.Res.XXXXX.h5
   wsa.gam.h5


The ``gamhelio.x.helio_serial.out`` contains the terminal output generated by the heliosphere version of the ``gamera`` executable (\ ``gamhelio.x``\ ). The file ``helio_serial.o13122765`` (the digits will be different) contains the PBS log for the job.

We can now analyze the results of the model run. A simple analysis can be run using the utility ``run_helio_serial_checks.py``\ :

.. code-block::

   #!shell
   $ run_helio_serial_checks.py -v --runid=wsa
   Computing volume-integrated magnetic pressure.
   Volume-integrated magnetic pressure (SUM(Pb*dV), code units):
   At start: 103622.17348083727
   At end: 90669.31440751288


Your values for the volume-integrated magnetic pressure should be very close to these values. This is actually a dummy statistic - it is not scientifically useful, but the code illustrates how to access and manipulate the results using the ``kaiju`` software.

Finally, generate a quick-look plot illustrating the model results. For this case, the quick-look plot shows several plots of different data generated by the model:

.. code-block::

   #!shell
   $ create_helio_serial_quicklook.py -v --runid=wsa


This script will create the file ``qkpichelio.png``. It should look like this:


.. image:: https://bitbucket.org/repo/kMoBzBp/images/2400191851-qkpichelio-min.png
   :target: https://bitbucket.org/repo/kMoBzBp/images/2400191851-qkpichelio-min.png
   :alt: qkpichelio.png


Running a test case (MPI)
^^^^^^^^^^^^^^^^^^^^^^^^^

3-D blast wave
~~~~~~~~~~~~~~

The first MPI-based example is the ``bw3d`` model. This is a simple 3-D model of a blast wave with an applied magnetic field. This test case is discussed at `Blast Wave Test <https://www.astro.princeton.edu/~jstone/Athena/tests/blast/blast.html>`_.

Begin by adding the directory containing the MPI kaiju binaries to your ``PATH`` environment variable.

.. code-block::

   #!shell
   export PATH=$KAIJUHOME/build_mpi/bin:$PATH


Next, add the directory containing the the ``bw3d`` example scripts to your ``PATH`` environment variable:

.. code-block::

   #!shell
   export PATH=$KAIJUHOME/quickstart/bw3d:$PATH


Create a new directory to run the ``bw3d`` example (it can be anywhere).

.. code-block::

   #!shell
   cd $HOME
   mkdir -p kaiju_test/bw3d
   cd kaiju_test/bw3d


The next step is to generate the configuration files for the ``bw3d`` model. This is done using the ``prepare_bw3d.py`` utility:

.. code-block::

   #!shell
   prepare_bw3d.py -v


This command will create 3 files in your current directory:


* ``bw3d.ini`` - An ``.ini``\ -format initialization file for running ``gamera_mpi.x`` on the ``bw3d`` example.
* ``bw3d.pbs`` - An PBS job script for running ``gamera_mpi.x`` on the ``bw3d`` example.
* ``bw3d.xml`` - An XML-format initialization file for running ``gamera_mpi.x`` on the ``bw3d`` example, created from ``bw3d.ini``.

NOTE: The conversion from ``.ini`` to ``.xml`` is still under development, so both initialization files are created from templates.

The PBS script is designed for use on the Pleiades system at NASA Ames. If you are working in a non-HPC environment, the commands listed in this file can be executed manually on your command line. On Pleiades, the job is submitted using the ``qsub`` command:

.. code-block::

   #!shell
   qsub bw3d.pbs


Once the job has been accepted in the queue, the run should take about an hour (on Pleiades).

When complete, you should see in your working directory the output files created by ``gamera_mpi.x``\ :

.. code-block::

   #!shell
   bash-4.2$ ls -1
   bw3d_0002_0002_0002_0000_0000_0000.gam.h5
   bw3d_0002_0002_0002_0000_0000_0000.gam.Res.00000.h5
   bw3d_0002_0002_0002_0000_0000_0000.gam.Res.XXXXX.h5
   bw3d_0002_0002_0002_0000_0000_0001.gam.h5
   bw3d_0002_0002_0002_0000_0000_0001.gam.Res.00000.h5
   bw3d_0002_0002_0002_0000_0000_0001.gam.Res.XXXXX.h5
   bw3d_0002_0002_0002_0000_0001_0000.gam.h5
   bw3d_0002_0002_0002_0000_0001_0000.gam.Res.00000.h5
   bw3d_0002_0002_0002_0000_0001_0000.gam.Res.XXXXX.h5
   bw3d_0002_0002_0002_0000_0001_0001.gam.h5
   bw3d_0002_0002_0002_0000_0001_0001.gam.Res.00000.h5
   bw3d_0002_0002_0002_0000_0001_0001.gam.Res.XXXXX.h5
   bw3d_0002_0002_0002_0001_0000_0000.gam.h5
   bw3d_0002_0002_0002_0001_0000_0000.gam.Res.00000.h5
   bw3d_0002_0002_0002_0001_0000_0000.gam.Res.XXXXX.h5
   bw3d_0002_0002_0002_0001_0000_0001.gam.h5
   bw3d_0002_0002_0002_0001_0000_0001.gam.Res.00000.h5
   bw3d_0002_0002_0002_0001_0000_0001.gam.Res.XXXXX.h5
   bw3d_0002_0002_0002_0001_0001_0000.gam.h5
   bw3d_0002_0002_0002_0001_0001_0000.gam.Res.00000.h5
   bw3d_0002_0002_0002_0001_0001_0000.gam.Res.XXXXX.h5
   bw3d_0002_0002_0002_0001_0001_0001.gam.h5
   bw3d_0002_0002_0002_0001_0001_0001.gam.Res.00000.h5
   bw3d_0002_0002_0002_0001_0001_0001.gam.Res.XXXXX.h5
   bw3d.ini
   bw3d.o12920469
   bw3d.pbs
   bw3d.xml
   gamera_mpi.x.bw3d.out


The ``gamera_mpi.x.bw3d.out`` contains the terminal output generated by the ``gamera_mpi.x`` executable. The file ``bw3d.o12920469`` (the digits will be different) contains the PBS log for the job.

We can now analyze the results of the model run. A simple analysis can be run using the utility ``run_bw3d_checks.py``\ :

.. code-block::

   #!shell
   $ run_bw3d_checks.py -v
   Asymmetry metric (SUM(ABS(Pb - Pb.T)*dV), code units):
   At start: 0.0
   At end: 0.001626661792365313


Your values for the asymmetry metric should be very close to these values.

Finally, generate a quick-look plot illustrating the model results. For this case, the quick-look plot shows the magnetic pressure in the first and last simulation frames:

.. code-block::

   #!shell
   $ create_bw3d_quicklook.py -v


This script will create the file ``bw3d_quicklook.png``. It should look like this:


.. image:: https://bitbucket.org/repo/kMoBzBp/images/431662986-bw3d_quicklook.png
   :target: https://bitbucket.org/repo/kMoBzBp/images/431662986-bw3d_quicklook.png
   :alt: bw3d_quicklook.png


Heliosphere example (MPI)
~~~~~~~~~~~~~~~~~~~~~~~~~

The next example is an MPI version of the serial heliosphere model described above. The mechanics of the MPI version are slightly different, but the results should be identical.

Begin by adding the directory containing the MPI kaiju binaries to your ``PATH`` environment variable.

.. code-block::

   #!shell
   export PATH=$KAIJUHOME/build_mpi/bin:$PATH


Next, add the directory containing the the ``helio_mpi`` example scripts to your ``PATH`` environment variable:

.. code-block::

   #!shell
   export PATH=$KAIJUHOME/quickstart/helio_mpi:$PATH


The next step is to generate the configuration and input files for the ``helio_mpi`` model. This is done using the ``prepare_helio_mpi.py`` utility:

.. code-block::

   #!shell
   prepare_helio_mpi.py -v


This command will create the following files in your current directory:


* ``helio_mpi.ini`` - An ``.ini-format`` initialization file for running ``gamera_mpi.x`` on the ``helio_mpi`` example.
* ``helio_mpi.pbs`` - A PBS job script for running ``gamera_mpi.x`` on the ``helio_mpi`` example.
* ``helio_mpi.xml`` - An XML-format initialization file for running ``gamera_mpi.x`` on the ``helio_mpi`` example, created from ``helio_mpi.ini``.
* heliogrid.h5 - A HDF5 file containing the grid to use for the simulation.
* innerbc.h5 - A HDF5 file containing the inner boundary conditions derived from the the WSA (Wang-Sheeley-Arge) model used for this example.

NOTE: The detailed steps for manually generating these files are described `here <helioQuickStart>`_.

NOTE: The conversion from ``.ini`` to ``.xml`` is still under development, so both initialization files are created from templates.

The PBS script is designed for use on the Pleiades system at NASA Ames. If you are working in a non-HPC environment, the commands listed in this file can be executed manually on your command line. On Pleiades, the job is submitted using the ``qsub`` command:

.. code-block::

   #!shell
   qsub helio_mpi.pbs


Once the job has been accepted in the queue, the run should take about 30 minutes (on Pleiades).

When complete, you should see in your working directory the output files created by ``gamera_mpi.x``\ , along with the files created by ``prepare_helio_mpi.py``\ :

.. code-block::

   #!shell
   bash-4.2$ ls -1
   C02DT5CZMD6T-ML:20220328 winteel1$ ls -1 helio_mpi
   gamera_mpi.x.helio_mpi.out
   helio_mpi.ini
   helio_mpi.o13124130
   helio_mpi.pbs
   helio_mpi.xml
   heliogrid.h5
   innerbc.h5
   qkpichelio.png
   wsa_0004_0002_0004_0000_0000_0000.gam.Res.00000.h5
   wsa_0004_0002_0004_0000_0000_0000.gam.Res.XXXXX.h5
   wsa_0004_0002_0004_0000_0000_0000.gam.h5
   wsa_0004_0002_0004_0000_0000_0001.gam.Res.00000.h5
   wsa_0004_0002_0004_0000_0000_0001.gam.Res.XXXXX.h5
   wsa_0004_0002_0004_0000_0000_0001.gam.h5
   wsa_0004_0002_0004_0000_0000_0002.gam.Res.00000.h5
   wsa_0004_0002_0004_0000_0000_0002.gam.Res.XXXXX.h5
   wsa_0004_0002_0004_0000_0000_0002.gam.h5
   wsa_0004_0002_0004_0000_0000_0003.gam.Res.00000.h5
   wsa_0004_0002_0004_0000_0000_0003.gam.Res.XXXXX.h5
   wsa_0004_0002_0004_0000_0000_0003.gam.h5
   wsa_0004_0002_0004_0000_0001_0000.gam.Res.00000.h5
   wsa_0004_0002_0004_0000_0001_0000.gam.Res.XXXXX.h5
   wsa_0004_0002_0004_0000_0001_0000.gam.h5
   wsa_0004_0002_0004_0000_0001_0001.gam.Res.00000.h5
   wsa_0004_0002_0004_0000_0001_0001.gam.Res.XXXXX.h5
   wsa_0004_0002_0004_0000_0001_0001.gam.h5
   wsa_0004_0002_0004_0000_0001_0002.gam.Res.00000.h5
   wsa_0004_0002_0004_0000_0001_0002.gam.Res.XXXXX.h5
   wsa_0004_0002_0004_0000_0001_0002.gam.h5
   wsa_0004_0002_0004_0000_0001_0003.gam.Res.00000.h5
   wsa_0004_0002_0004_0000_0001_0003.gam.Res.XXXXX.h5
   wsa_0004_0002_0004_0000_0001_0003.gam.h5
   wsa_0004_0002_0004_0001_0000_0000.gam.Res.00000.h5
   wsa_0004_0002_0004_0001_0000_0000.gam.Res.XXXXX.h5
   wsa_0004_0002_0004_0001_0000_0000.gam.h5
   wsa_0004_0002_0004_0001_0000_0001.gam.Res.00000.h5
   wsa_0004_0002_0004_0001_0000_0001.gam.Res.XXXXX.h5
   wsa_0004_0002_0004_0001_0000_0001.gam.h5
   wsa_0004_0002_0004_0001_0000_0002.gam.Res.00000.h5
   wsa_0004_0002_0004_0001_0000_0002.gam.Res.XXXXX.h5
   wsa_0004_0002_0004_0001_0000_0002.gam.h5
   wsa_0004_0002_0004_0001_0000_0003.gam.Res.00000.h5
   wsa_0004_0002_0004_0001_0000_0003.gam.Res.XXXXX.h5
   wsa_0004_0002_0004_0001_0000_0003.gam.h5
   wsa_0004_0002_0004_0001_0001_0000.gam.Res.00000.h5
   wsa_0004_0002_0004_0001_0001_0000.gam.Res.XXXXX.h5
   wsa_0004_0002_0004_0001_0001_0000.gam.h5
   wsa_0004_0002_0004_0001_0001_0001.gam.Res.00000.h5
   wsa_0004_0002_0004_0001_0001_0001.gam.Res.XXXXX.h5
   wsa_0004_0002_0004_0001_0001_0001.gam.h5
   wsa_0004_0002_0004_0001_0001_0002.gam.Res.00000.h5
   wsa_0004_0002_0004_0001_0001_0002.gam.Res.XXXXX.h5
   wsa_0004_0002_0004_0001_0001_0002.gam.h5
   wsa_0004_0002_0004_0001_0001_0003.gam.Res.00000.h5
   wsa_0004_0002_0004_0001_0001_0003.gam.Res.XXXXX.h5
   wsa_0004_0002_0004_0001_0001_0003.gam.h5
   wsa_0004_0002_0004_0002_0000_0000.gam.Res.00000.h5
   wsa_0004_0002_0004_0002_0000_0000.gam.Res.XXXXX.h5
   wsa_0004_0002_0004_0002_0000_0000.gam.h5
   wsa_0004_0002_0004_0002_0000_0001.gam.Res.00000.h5
   wsa_0004_0002_0004_0002_0000_0001.gam.Res.XXXXX.h5
   wsa_0004_0002_0004_0002_0000_0001.gam.h5
   wsa_0004_0002_0004_0002_0000_0002.gam.Res.00000.h5
   wsa_0004_0002_0004_0002_0000_0002.gam.Res.XXXXX.h5
   wsa_0004_0002_0004_0002_0000_0002.gam.h5
   wsa_0004_0002_0004_0002_0000_0003.gam.Res.00000.h5
   wsa_0004_0002_0004_0002_0000_0003.gam.Res.XXXXX.h5
   wsa_0004_0002_0004_0002_0000_0003.gam.h5
   wsa_0004_0002_0004_0002_0001_0000.gam.Res.00000.h5
   wsa_0004_0002_0004_0002_0001_0000.gam.Res.XXXXX.h5
   wsa_0004_0002_0004_0002_0001_0000.gam.h5
   wsa_0004_0002_0004_0002_0001_0001.gam.Res.00000.h5
   wsa_0004_0002_0004_0002_0001_0001.gam.Res.XXXXX.h5
   wsa_0004_0002_0004_0002_0001_0001.gam.h5
   wsa_0004_0002_0004_0002_0001_0002.gam.Res.00000.h5
   wsa_0004_0002_0004_0002_0001_0002.gam.Res.XXXXX.h5
   wsa_0004_0002_0004_0002_0001_0002.gam.h5
   wsa_0004_0002_0004_0002_0001_0003.gam.Res.00000.h5
   wsa_0004_0002_0004_0002_0001_0003.gam.Res.XXXXX.h5
   wsa_0004_0002_0004_0002_0001_0003.gam.h5
   wsa_0004_0002_0004_0003_0000_0000.gam.Res.00000.h5
   wsa_0004_0002_0004_0003_0000_0000.gam.Res.XXXXX.h5
   wsa_0004_0002_0004_0003_0000_0000.gam.h5
   wsa_0004_0002_0004_0003_0000_0001.gam.Res.00000.h5
   wsa_0004_0002_0004_0003_0000_0001.gam.Res.XXXXX.h5
   wsa_0004_0002_0004_0003_0000_0001.gam.h5
   wsa_0004_0002_0004_0003_0000_0002.gam.Res.00000.h5
   wsa_0004_0002_0004_0003_0000_0002.gam.Res.XXXXX.h5
   wsa_0004_0002_0004_0003_0000_0002.gam.h5
   wsa_0004_0002_0004_0003_0000_0003.gam.Res.00000.h5
   wsa_0004_0002_0004_0003_0000_0003.gam.Res.XXXXX.h5
   wsa_0004_0002_0004_0003_0000_0003.gam.h5
   wsa_0004_0002_0004_0003_0001_0000.gam.Res.00000.h5
   wsa_0004_0002_0004_0003_0001_0000.gam.Res.XXXXX.h5
   wsa_0004_0002_0004_0003_0001_0000.gam.h5
   wsa_0004_0002_0004_0003_0001_0001.gam.Res.00000.h5
   wsa_0004_0002_0004_0003_0001_0001.gam.Res.XXXXX.h5
   wsa_0004_0002_0004_0003_0001_0001.gam.h5
   wsa_0004_0002_0004_0003_0001_0002.gam.Res.00000.h5
   wsa_0004_0002_0004_0003_0001_0002.gam.Res.XXXXX.h5
   wsa_0004_0002_0004_0003_0001_0002.gam.h5
   wsa_0004_0002_0004_0003_0001_0003.gam.Res.00000.h5
   wsa_0004_0002_0004_0003_0001_0003.gam.Res.XXXXX.h5
   wsa_0004_0002_0004_0003_0001_0003.gam.h5


The ``gamera_mpi.x.helio_mpi.out`` contains the terminal output generated by the heliosphere version of the ``gamera_mpi.x`` executable. The file ``helio_mpi.o13124130`` (the digits will be different) contains the PBS log for the job.

We can now analyze the results of the model run. A simple analysis can be run using the utility ``run_helio_mpi_checks.py``\ :

.. code-block::

   #!shell
   $ run_helio_mpi_checks.py -v --runid=wsa
   Computing volume-integrated magnetic pressure.
   Volume-integrated magnetic pressure (SUM(Pb*dV), code units):
   At start: 103622.17348083727
   At end: 90669.31440751288


Your values for the volume-integrated magnetic pressure should be very close to these values. This is actually a dummy statistic - it is not scientifically useful, but the code illustrates how to access and manipulate the results using the ``kaiju`` software.

Finally, generate a quick-look plot illustrating the model results. For this case, the quick-look plot shows several plots of different data generated by the model:

.. code-block::

   #!shell
   $ create_helio_mpi_quicklook.py -v --runid=wsa


This script will create the file ``qkpichelio.png``. It should look like this:


.. image:: https://bitbucket.org/repo/kMoBzBp/images/4029922679-qkpichelio-min.png
   :target: https://bitbucket.org/repo/kMoBzBp/images/4029922679-qkpichelio-min.png
   :alt: qkpichelio-min.png

