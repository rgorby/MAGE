Wave model input for RCM diffuse precipitation
=====================

To run MAGE with wave model input for the diffuse precipitation:

Make sure the parameters for the polynomial fit, chorus_polynomial.txt, is inside kaipy/rcm/wmutils

To generate rcmconfig.h5 that contains the electron lifetime based on the wave model,
just run kaiju/scripts/preproc/genRCM.py. The wave model is added by default.

To generate rcmconfig.h5 without the wave model

.. code-block:: python

   python genRCM.py --noWaveModel

To add the electron lifetime to an existing rcmconfig.h5 file,
run kaiju/scripts/preproc/genRCM.py with 'addWM' option on and enter input file name.

.. code-block:: python

   python genRCM.py --addWM True -i rcmconfig.h5

To set the maximum Kp index allowed in the chorus wave model (maxKp <= 6, default 6)

.. code-block:: python

   python genRCM.py -maxKp 6

The generated rcmconfig.h5 should contain the arrays for the wave model as follows:

* Eki                      Dataset {155}
* Kpi                      Dataset {6}
* Li                       Dataset {41}
* MLTi                     Dataset {97}
* Taui                     Dataset {155, 41, 97, 6}

In the XML file, to run MAGE with wave model input (Dedong Wang Chorus + Orlova16 Hiss), set the loss method in the RCM section as

.. code-block:: python

   <RCM>
       <loss eLossMethod="WM"/>
   </RCM>

Go to Step 1. if you see“Wave model is missing in rcmconfig.h5” in .out file of the MAGE run.
