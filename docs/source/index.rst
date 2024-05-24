.. role:: raw-html-m2r(raw)
   :format: html


Welcome to the kaiju documentation
==================================

Use of any code contained in this repository or data produced by it is expected
to respect these roadrules_.

..  For questions, please, `contact us <mailto:MAGEEC@LISTSERV.JHUAPL.EDU>`_.

.. **Code improvements and contributions are welcome, and should follow the
.. pull request process outlined on `this page <userGuide/contributingGuide>`_.**

.. **If you encounter issues go through the 
.. `Troubleshooting Guide <quickStart/trouble>`_ and `FAQ <userGuide/FAQ>`_ as well as `Google <https://www.google.com/>`_ before seeking help about your
.. issue. If you're still having difficulty post a detailed description of
.. your problem in the `#kaijuhelp <https://nasa-drive-cgs.slack.com/archives/C011V6V7YSJ>`_ Slack channel on the CGS workspace,
.. where the developers and other users can try to assist you.**

.. ----

.. **Quick Start Guide** :raw-html-m2r:`<a name="quickStart"></a>`
.. =======================================================================

.. Building the code
.. -----------------


.. * `Prerequisites <quickStart/prerequisites>`_ - Setting up your system before you
..   clone the git repo
.. * `Building kaiju <building/build>`_ - Building the code (new version)
.. * `Verify installation <quickStart/tests>`_ - Testing your installation
.. * `Performance <userGuide/performance>`_ - Performance Considerations

.. Running the code
.. ----------------


.. * `Geospace Quick Start Guide <quickStart/geoQuickStart>`_ - Starting guide for
..   geospace simulations
.. * `Heliosphere Quick Start Guide <quickStart/helioQuickStart>`_ - Starting guide for
..   heliospheric simulations
.. * `CHIMP Quick Start <quickStart/chimpQuickStart>`_ - Starting guide for test particle
..   simulations
.. * `Planetary Magnetosphere Quick Start <quickStart/planetaryQuickStart>`_ - Starting
..   guide for non-terrestrial magnetospheric simulations

.. Viewing the results
.. -------------------


.. * `Quick look plots <quickStart/quickLook>`_

.. Distributing the ``kaiju`` Python code (\ ``kaipy``\ ):
.. -----------------------------------------------------------


.. * `Building a pip-installable package <distributing/kaipy_for_pip>`_
.. * `Building a conda-installable package <distributing/kaipy_for_conda>`_

.. ----

.. **User Guide** :raw-html-m2r:`<a name="userGuide"></a>`
.. ===============================================================

.. User Rules
.. ----------


.. * `Contributing Guide <userGuide/contributingGuide>`_
.. * `Wiki Contributing Guide <userGuide/wikiContributing>`_
.. * `Development Roadmap <userGuide/developmentRoadmap>`_
.. * `FAQ <userGuide/FAQ>`_

.. Magnetosphere Simulations with MAGE
.. -----------------------------------


.. * `GAMERA REMIX <userGuide/gameraRemix>`_ - Running a GAMERA REMIX simulation
.. * `GAMERA RCM <userGuide/gameraRCM>`_ - Running coupled GAMERA RCM
.. * `MAGE <userGuide/mage>`_ - Running MAGE 1.0 (GAMERA + RCM + TIEGCM)
.. * `HIDRA <userGuide/hidra>`_ - Running Hidra ionospheric outflow simulations
.. * XML Files - configuration, defaults, and description

..   * `Generating XML <userGuide/generatingXML>`_ - Useful info on generating XML files
..   * `GAMERA <userGuide/gameraXML>`_ - Basic GAMERA info
..   * `VOLTRON <userGuide/voltronXML>`_ - VOLTRON info

.. * MAGE-FAQ - Frequently asked questions
.. * `Computational Costs <userGuide/compCosts>`_ - Costs for running various systems
.. * Analysis tools

..   * `Ground Magnetic Field Calculations <userGuide/groundMag>`_ - Biot-Savart calculation of magnetic field perturbations
..   * `SuperMage <userGuide/superMAGE>`_ - Data Model Comparison with SuperMag

.. * `Outer planets and exoplanets <userGuide/exoOuterPlanets>`_ - Non-Earth planetary magnetospheres

.. Heliosphere Simulations with GAMERA-Helio
.. -----------------------------------------


.. * `Steady State Run <userGuide/steadyStateRun>`_ - Creating grid, boundary conditions and submit a mpi gamera helio run
.. * `Helio CME <userGuide/helio-cme>`_ - Running GAMERA-Helio with an embedded CME
.. * `Gibson-Low CME Model <userGuide/gibson-low>`_ - Running standalone Gibson-Low CME Model
.. * XML Files - configuration, defaults, and description
.. * `Visualizing the results <userGuide/helioVisualizing>`_
.. * `Computational Costs <userGuide/compCostsHelio>`_ - Costs for running Helio simulations

.. Particle Simulations and Analysis Tools with CHIMP
.. --------------------------------------------------


.. * `Test particle integration: ``push.x`` <userGuide/push.x>`_ - explanation of XML parameters along with an example XML file and job script
.. * `Phase space density calculation: ``psd.x`` <userGuide/psd.x>`_
.. * `Extract 2D slice: ``slice.x`` <userGuide/slice.x>`_
.. * `Extract 3D  subdomain: ``chop.x`` <userGuide/chop.x>`_
.. * `Magnetic field line tracer: ``trace.x`` <userGuide/trace.x>`_
.. * `Ground magnetic perturbation calculation calculation: ``calcdb.x`` <userGuide/groundMag>`_
.. * `CHIMP XML <userGuide/chimpXML>`_ - explanation of XML input files for CHIMP

.. Rice Convection Model
.. ---------------------


.. * `RCMX <userGuide/rcmx>`_ - Stand-alone RCM driver (rcm.x)
.. * `XML Files <userGuide/rcm_xml>`_ - configuration, defaults, and description
.. * `lambdautils <userGuide/lambdaUtils>`_ - Generating and testing rcmconfigs with lambdautils
.. * `Diffuse precipitation <userGuide/wmutils>`_ - Electron loss and precipitation prescribed by wave models

.. Comparing results with satellite data
.. -------------------------------------


.. * `msphsatcomp <userGuide/msphsatcomp>`_ - Comparison of terrestrial magnetosphere runs with satellite data
.. * `heliosatcomp <userGuide/heliosatcomp>`_ - Comparison of heliosphere runs with satellite data

.. Tools
.. -----


.. * `Containerization <Pokeball>`_ - Early work on containerization/cloud computing
.. * `Sublime Text <userGuide/Sublime_Text>`_ - Using sublime text for development
.. * `Debugging <userGuide/debugging>`_ - Setting up and running the Allinea parallel debugger
.. * `Globus Endpoint <userGuide/globusEndpoint>`_ - Sharing data at NCAR via Globus Endpoint
.. * `NASA JupyterLab <userGuide/nasaJupyter>`_ - Using JupyterLab on NASA HEC systems to analyze results

.. Code information
.. ----------------


.. * `Repository Organization <userGuide/codeOrg>`_  -  Structure of the kaiju repository
.. * `Interpolation function in MAGE <userGuide/interpolationInMAGE>`_ - Info on how interpolation is done in MAGE
.. * `Derivation of precipitation <userGuide/derivation_of_precipitation>`_ - Derivation of auroral precipitation in MAGE.

.. Testing
.. -------


.. * `Unit Testing <userGuide/unitTesting>`_ - Setting up unit testing framework for Kaiju codes
.. * `Adding new tests <userGuide/addingNewTests>`_ - Steps for adding new tests
.. * `TIE-GCM Benchmarks <userGuide/tiegcmBenchmarks>`_ - Instructions for running the TIEGCM benchmarks
