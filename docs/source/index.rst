Welcome to the kaiju documentation
==================================

**Read this First!**
--------------------

Use of any code contained in this repository or data produced by it is expected to respect these :doc:`rules of the road <roadrules>`. For questions, please, `contact us <MAGEEC@LISTSERV.JHUAPL.EDU>`_.

Code improvements and contributions are welcome, and should follow the pull request process outlined on :doc:`this page <userGuide/contributingGuide>`.

If you encounter issues go through the :doc:`Troubleshooting Guide <quickStart/trouble>` and :doc:`FAQ <userGuide/FAQ>` as well as `Google <https://www.google.com/>`_ before seeking help about your issue. If you're still having difficulty post a detailed description of your problem in the `#kaijuhelp <https://nasa-drive-cgs.slack.com/archives/C011V6V7YSJ>`_ Slack channel on the CGS workspace, where the developers and other users can try to assist you.

----

**Quick Start Guide**
---------------------

Building the code
~~~~~~~~~~~~~~~~~
* :doc:`Prerequisites <quickStart/prerequisites>` - Setting up your system before you clone the git repo
* :doc:`Building kaiju <building/build>` - Building the code (new version)
* :doc:`Verify installation <quickStart/tests>` - Testing your installation
* :doc:`Performance <userGuide/performance>` - Performance Considerations

Running the code
~~~~~~~~~~~~~~~~
* :doc:`Geospace Quick Start Guide <quickStart/geoQuickStart>` - Starting guide for geospace simulations
* :doc:`Heliosphere Quick Start Guide <quickStart/helioQuickStart>` - Starting guide for heliospheric simulations
* :doc:`CHIMP Quick Start <quickStart/chimpQuickStart>` - Starting guide for test particle simulations
* :doc:`Planetary Magnetosphere Quick Start <quickStart/planetaryQuickStart>` - Starting guide for non-terrestrial magnetospheric simulations

Viewing the results
~~~~~~~~~~~~~~~~~~~
* :doc:`Quick look plots <quickStart/quickLook>`

Distributing the ``kaiju`` Python code (``kaipy``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* :doc:`Building a pip-installable package <distributing/kaipy_for_pip>`
* :doc:`Building a conda-installable package <distributing/kaipy_for_conda>`

----

**User Guide**
--------------

User Rules
~~~~~~~~~~
* :doc:`Contributing Guide <userGuide/contributingGuide>`
* :doc:`Wiki Contributing Guide <userGuide/wikiContributing>`
* :doc:`Development Roadmap <userGuide/developmentRoadmap>`
* :doc:`FAQ <userGuide/FAQ>`

Magnetosphere Simulations with MAGE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* :doc:`GAMERA REMIX <userGuide/gameraRemix>` - Running a GAMERA REMIX simulation
* :doc:`GAMERA RCM <userGuide/gameraRCM>` - Running coupled GAMERA RCM
* :doc:`MAGE <userGuide/mage>` - Running MAGE 1.0 (GAMERA + RCM + TIEGCM)
* :doc:`HIDRA <userGuide/hidra>` - Running Hidra ionospheric outflow simulations
* XML Files - configuration, defaults, and description
    * :doc:`Generating XML <userGuide/generatingXML>` - Useful info on generating XML files
    * :doc:`GAMERA <userGuide/gameraXML>` - Basic GAMERA info
    * :doc:`VOLTRON <userGuide/voltronXML>` - VOLTRON info
* MAGE-FAQ - Frequently asked questions
* :doc:`Computational Costs <userGuide/compCosts>` - Costs for running various systems
* Analysis tools
    * :doc:`Ground Magnetic Field Calculations <userGuide/groundMag>` - Biot-Savart calculation of magnetic field perturbations
    * :doc:`SuperMage <userGuide/superMAGE>` - Data Model Comparison with SuperMag
* :doc:`Outer planets and exoplanets <userGuide/exoOuterPlanets>` - Non-Earth planetary magnetospheres

Heliosphere Simulations with GAMERA-Helio
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* :doc:`Steady State Run <userGuide/steadyStateRun>` - Creating grid, boundary conditions and submit a mpi gamera helio run
* :doc:`Helio CME <userGuide/helio-cme>` - Running GAMERA-Helio with an embedded CME
* :doc:`Gibson-Low CME Model <userGuide/gibson-low>` - Running standalone Gibson-Low CME Model
* XML Files - configuration, defaults, and description
* :doc:`Visualizing the results <userGuide/helioVisualizing>`
* :doc:`Computational Costs <userGuide/compCostsHelio>` - Costs for running Helio simulations

Particle Simulations and Analysis Tools with CHIMP
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* :doc:`Test particle integration: push.x <userGuide/push.x>` - explanation of XML parameters along with an example XML file and job script
* :doc:`Phase space density calculation: psd.x <userGuide/psd.x>`
* :doc:`Extract 2D slice: slice.x <userGuide/slice.x>`
* :doc:`Extract 3D  subdomain: chop.x <userGuide/chop.x>`
* :doc:`Magnetic field line tracer: trace.x <userGuide/trace.x>`
* :doc:`Ground magnetic perturbation calculation calculation: calcdb.x <userGuide/groundMag>`
* :doc:`CHIMP XML <userGuide/chimpXML>` - explanation of XML input files for CHIMP

Rice Convection Model
~~~~~~~~~~~~~~~~~~~~~
* :doc:`RCMX <userGuide/rcmx>` - Stand-alone RCM driver (rcm.x)
* :doc:`XML Files <userGuide/rcm_xml>` - configuration, defaults, and description
* :doc:`lambdautils <userGuide/lambdaUtils>` - Generating and testing rcmconfigs with lambdautils
* :doc:`Diffuse precipitation <userGuide/wmutils>` - Electron loss and precipitation prescribed by wave models

Comparing results with satellite data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* :doc:`msphsatcomp <userGuide/msphsatcomp>` - Comparison of terrestrial magnetosphere runs with satellite data
* :doc:`heliosatcomp <userGuide/heliosatcomp>` - Comparison of heliosphere runs with satellite data

Tools
~~~~~
* :doc:`Containerization <Pokeball>` - Early work on containerization/cloud computing
* :doc:`Sublime Text <userGuide/Sublime_Text>` - Using sublime text for development
* :doc:`Debugging <userGuide/debugging>` - Setting up and running the Allinea parallel debugger
* :doc:`Globus Endpoint <userGuide/globusEndpoint>` - Sharing data at NCAR via Globus Endpoint
* :doc:`NASA JupyterLab <userGuide/nasaJupyter>` - Using JupyterLab on NASA HEC systems to analyze results

Code information
~~~~~~~~~~~~~~~~
* :doc:`Repository Organization <userGuide/codeOrg>`  -  Structure of the kaiju repository
* :doc:`Interpolation function in MAGE <userGuide/interpolationInMAGE>` - Info on how interpolation is done in MAGE
* :doc:`Derivation of precipitation <userGuide/derivation_of_precipitation>` - Derivation of auroral precipitation in MAGE.

Testing
~~~~~~~
* :doc:`Unit Testing <userGuide/unitTesting>` - Setting up unit testing framework for Kaiju codes
* :doc:`Adding new tests <userGuide/addingNewTests>` - Steps for adding new tests
* :doc:`TIE-GCM Benchmarks <userGuide/tiegcmBenchmarks>` - Instructions for running the TIEGCM benchmarks
