Comparison of heliosphere simulation results with satellite data
================================================================

Introduction
------------

Comparison of simulation results to observed data is a critical part of model
validation. The MAGE team has developed a set of tools to facilitate
comparison of the results of heliospheric simulations with measured spacecraft
data. The initial version of these tools supports comparison of the results
from runs of ``gamhelio.x`` or ``gamhelio_mpi.x`` to data from the
`ACE spacecraft <https://solarsystem.nasa.gov/missions/ace/in-depth/>`_. The
observations for the desired time period are retrieved from
`CDAWeb <https://cdaweb.gsfc.nasa.gov/>`_. Other spacecraft and data sources
will be added in the future as time permits.

The satellite comparison script
-------------------------------

The basic tool for heliospheric comparisons is the script ``helioSatComp.py``.
This script is available under the ``kaiju/scripts/datamodel`` directory of
your clone of the ``kaiju`` repository. Executing this script with the ``-h``
option will provide the following help text:

.. code-block:: bash

   usage: helioSatComp.py [-h] [-c command] [-d] [--deltaT deltaT] [-i runid]
                          [-k] [-n number_segments] [-p path] [--Rin Rin]
                          [--Rout Rout] [-s satellite_id] [-v]

   Extract satellite trajectory and observations for various
   heliospheric spacecraft from CDAWeb. Produce comparisons between the
   observations and corresponding gamhelio model results.

   optional arguments:
     -h, --help            show this help message and exit
     -c command, --cmd command
                           Full path to sctrack.x command (default: $KAIJUHOME/build/bin/sctrack.x).
     -d, --debug           Print debugging output (default: False).
     --deltaT deltaT       Time interval (seconds) for ephemeris points from CDAWeb (default: 3600.0).
     -i runid, --id runid  ID string of the run (default: hsphere)
     -k, --keep            Keep intermediate files (default: False).
     -n number_segments, --numSeg number_segments
                           Number of segments to simultaneously process (default: 1).
     -p path, --path path  Path to directory containing gamhelio results (default: .)
     --Rin Rin             Inner boundary (Rsun) of gamhelio grid (default: 21.5).
     --Rout Rout           Outer boundary (Rsun) of gamhelio grid (default: 220.0).
     -s satellite_id, --satId satellite_id
                           Name of Satellite to compare
     -v, --verbose         Print verbose output (default: False).


The ``-h, --help`` argument will print the help message above.

The ``-c command, --cmd command`` argument specifies the path to the Fortran
program ``sctrack.x``, which is used to interpolate the simulation results
to the times and locations for the ephemeris points for the specified
satellite. If you have added ``$KAIJUHOME/build/bin`` to your ``PATH``
environment variable, you should not have to set this option.

The ``-d, --debug`` argument causes the script to output a variety of
debugging output.

The ``--deltaT deltaT`` argument allows you to specify the time interval for
the satellite data that will be retrieved from CDAWeb. The default value of
3600 seconds (1 hour) is usually sufficient, but smaller values of this
argument will allow more detailed comparisons.

The ``-i runid, --id runid`` argument allows you to specify the run ID of the
simulation results file you are using. The run ID is the identifying substring
of the file name. For example, in the simulation results file
``msphere.gam.h5``, the run ID is ``msphere``. Set this argument as needed for
the file you are using for your comparison.

The ``-k, --keep`` argument forces the script to retain intermediate files
generated during the data retrieval from CDAWeb. This is useful for debugging.

The ``-n number_segments, --numSeg number_segments`` allows you to specify the
number of segments in your simulation results. A segment is a chunk of the
ephemeris that is processed in parallel with other chunks. For example, if you
use the argument ``-n 2``, the ephemeris will be split into 2 separate
portions and processed in parallel.

The ``-p path, --path path`` argument allows you to specify the path to a
directory containing the simulation results file you wish to use for the
satellite comparison. By default, the current directory is searched for the
file.

The ``--Rin Rin`` and ``--Rout Rout`` arguments allow you to specify the inner
and outer radii (in solar radii) of your simulation results. You should not
need to change these values.

The ``-s satellite_id, --satId satellite_id`` argument allows you to specify
the name of the spacecraft to use for the data comparison. Currently, only
'ACE' is supported.

The ``-v, --verbose`` argument causes the script to print informational
messages during execution.

Procedure
---------

The basic procedure for the satellite comparison is simple:


Run the ``helioSatComp.py`` script.

Examine the resulting plots.

You will almost always have to specify at least a few arguments to the
``helioSatComp.py`` script in order to perform the comparison with satellite
data. This is best explained using an example.

Example
-------

For this example, assume we have a directory containing a set of 1 or more HDF
5 files generated by ``gamhelio.x`` or ``gamhelio_mpi.x``. It might look like
this, if generated by ``gamhelio.x``:

.. code-block:: bash

   heliogrid.h5    wsa.gam.Res.00000.h5.    wsa.gam.h5
   innerbc.h5      wsa.gam.Res.XXXXX.h5


The simulation results are in the file ``wsa.gam.h5``. To compare the
simulation results in this file to data measured by the ACE spacecraft during
the same time period, you will need to specify the run ID (``wsa``) and the
satellite ID (``ACE``). And for the sake of illustration, assume you have
built your code in a non-standard location (``$KAIJUHOME/build_serial`` rather
than ``$KAIJUHOME/build``). You can run the comparison with the following
command:

.. code-block:: bash

   helioSatComp.py -v --id=wsa --cmd=$KAIJUHOME/build_serial/bin/sctrack.x --satId=ACE


When this command completes, your directory will contain several new files:

.. code-block:: bash

   ACE-traj.png
   ACE-error.txt
   ACE.png
   ACE.comp.cdf
   sctrack.out


The file ``ACE-traj.png`` contains a set of plots illustrating the portion of
the spacecraft trajectory used in the comparison. An example is shown below:


.. image:: https://bitbucket.org/repo/kMoBzBp/images/3648821514-ACE-traj.png
   :target: https://bitbucket.org/repo/kMoBzBp/images/3648821514-ACE-traj.png
   :alt: ACE-traj.png

The file ``ACE-error.txt`` contains the error statistics for the variables
examined in the comparison. For the ACE spacecraft, these variables are:
``MagneticField``, ``Density``, ``Speed``, and ``Temperature``. Other
variables may be added later. This file looks something like this:

.. code-block:: bash

   Errors for: Density
   MAE: 12.838838182545178
   MSE: 171.86416943415404
   RMSE: 13.109697534045324
   MAPE: 0.9660806312821668
   RSE: 24.439823909708686
   PE: -23.439823909708686
   Errors for: Temperature
   MAE: 16256587733.017332
   MSE: 2.642806561220219e+20
   RMSE: 16256711110.246805
   MAPE: 214182.88316333634
   RSE: 44568596463.807556
   PE: -44568596462.807556
   Errors for: Speed
   MAE: 499.42446903122
   MSE: 260272.54722637442
   RMSE: 510.1691359013934
   MAPE: 0.9986149471382262
   RSE: 23.99323569578531
   PE: -22.99323569578531
   Errors for: Br
   MAE: 3.0560021664580166
   MSE: 12.410819281122311
   RMSE: 3.5228992720658807
   MAPE: 0.9957802458984243
   RSE: 1.0323378185670626
   PE: -0.03233781856706264

The file ``ACE.png`` contains plots which show the comparison of the
simulation and measured data. An example of this file is shown below.

.. image:: https://bitbucket.org/repo/kMoBzBp/images/1302341116-ACE.png
   :target: https://bitbucket.org/repo/kMoBzBp/images/1302341116-ACE.png
   :alt: ACE.png

The file ``ACE.comp.cdf`` is a `CDF file <https://cdf.gsfc.nasa.gov/>`_
containing the measured and interpolated simulated data for comparison.

The file ``sctrack.out`` contains a log generated by the program
``sctrack.x``, which performs the interpolation of the simulation results to
the ephemeris points.
