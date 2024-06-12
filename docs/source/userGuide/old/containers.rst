
Documentation related to containerizing and working with containers
===================================================================

TL;DR
=====

ContainerImages allow the user to pull a pre-built executable from a ContainerRegistry including all required softwares . Containers themselves are isolated processes controlled by a ContainerRuntime

PreRequisites to working with this containers
=============================================


* highlevel container runtime installed , e.g. ``docker``
* linux or Mac
* (easier on x86)

Supported UseCases
==================


* OpenSource: 


#. 

   * quickstart oneliners that allow a user to start working locally within (orders of) minutes

#. 

   * distribution of a controlled and versioned artefact for common architectures 


* OnBoarding: 


#. 

   * above plus:

#. 

   * development inside containers instead of directly on a laptop, to save debugging time

#. 

   * interactive server hosted environments that allow tutorials and live-session sharing

Getting started with containers - in general
============================================


* post some links here
* planned: self-service tutorial "containers for scientists/physicists"?

Supported Features
==================


* Paraview and Analysis: a container that includes all kaiju tooling to allow data analysis and rendering
* Local Development: a container that includes all of kaiju, such that user can run on (order of) 10 threads and modify the code/python without installing dependencies **
* Globus Authentication: ability to have a predefined function (or equivalent) that allows sign-on to Globus and access to Drag-and-Drop UI

** not open sourced (nor documented yet), code base currently very stale

Planned Features:
-----------------


* Secure remote mounts to avoid data downloads
* Merge code base (upgrade all packages)
* Browserbased rendering (e.g. WASM) and other fanciness

Current TODO (@Mike @Kareem @Brent @Constanze):
===============================================


* determine perfect time slot for merging code
* improve Globus Transfers/Authn, explore remote mounting via Globus
* data tiering on GCP (@C)
