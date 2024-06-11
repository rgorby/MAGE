
Overview and Quick Starts
=========================

Throughout the descriptions $KAIJUHOME refers to the base directory of the kaiju repository.


* [[kaitour|Kaiju Tour]] - Code organization and main data structures
* [[Quick Start]] - Simple guide for installing, compiling, and running a simple MHD test case
* [[Building on Pleiades]] - Information specific to building and running on NASA Pleiades
* [[Computational Costs]] - Tables with estimates for run times at various resolutions
* [[Gamerasphere]] - Quick start for setting up a terrestrial magnetosphere run
* [[CHIMP Quick Start]] - Starting guide for test particle simulations
* [[Generating XML Files]] - How to generate XML files for any setup, manually or with a helper python script
* [[Troubleshooting]] - Common problems and solutions
* [[FAQ]] - Frequently asked questions

----

Magnetosphere simulations w/ MAGE
=================================


* [[GAMERA-RCM]] - Running coupled Gamera+RCM
* [[VOLTRON XML]] - XML configuration, defaults, and description
* [[MAGE-FAQ]] - Frequently asked questions
* `Ground Magnetic Field Calculations <https://bitbucket.org/aplkaiju/kaiju/wiki/Ground%20Magnetometer%20Calculations>`_ - Biot-Savart calculation of magnetic field perturbations
* `SuperMage <https://bitbucket.org/aplkaiju/kaiju/wiki/Data%20Model%20Comparison%20with%20SuperMag>`_ - Data Model Comparison with SuperMag

----

GAMERA: MHD
-----------

Links to pages on running Gamera simulations and doing post-processing/visualization.


* [[GamXML]] - Explanation of Gamera XML input parameters
* [[GAMERA-FAQ]] - FAQ

GAMERA-MPI
~~~~~~~~~~

Links to pages on how to build and run MPI distributed Gamera simulations. These instructions should all be considered addendums to the detailed instructions provided above for these types of cases. Any case that can be run with MPI can also be run without it. Using MPI only requires a few additional XML parameters and changes to the job submission script.


* [[Building Gamera with MPI]]
* [[Gamera MPI]] - Basic MHD runs

GAMERA-Helio
------------

Links to pages on running Gamera for inner heliosphere simulations: grid, boundary conditions, post-processing.


* [[Steady state run]] - Creating grid, boundary conditions and submit a mpi gamera helio run
* [[Time-dependent run]] - Creating grid, boundary conditions and submit a mpi gamera helio run
* [[CME in the inner heliosphere]] - Coupling non-mpi Gamera-Helio and Gibson&Low CME model

CHIMP: Test-particles, field-line tracing, grid slicing, etc..
--------------------------------------------------------------

Links to pages on running CHIMP simulations and doing post-processing/visualization.


* [[CHIMP Quick Start]] - Starting guide for test particle simulations
* [[ChimpXML]] - Explanation of CHIMP XML input parameters

RCM: Rice Convection Model
--------------------------


* [[RCMX]] - Stand-alone RCM driver (rcm.x)
* [[lambdautils]] - Generating and testing rcmconfigs with lambdautils

----

Testing and Development
=======================


* [[Sublime Text]] - Using sublime text for development
* [[Unit Testing]] - Setting up unit testing framework for Kaiju codes.
* [[Adding New Tests]] - Adding new unit tests
* [[Development Roadmap]] - Next development steps
* [[Debugging]] - Setting up and running the Allinea parallel debugger

Wiki
====

Link to some wiki instructions [[Wiki Info]]
