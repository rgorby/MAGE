Frequently Asked Questions
==========================

.. This page contains frequently asked questions about running Kaiju.

.. Q1. I got ipo warnings #11021 of unresolved targets referenced in dynamic libraries when linking fortran executable gamera.x. The recipe failed. Example error message when "make gamera":
.. 'ipo: warning #11021: unresolved files_mp\ *checkandkill*
..         Referenced in libgamlib.a(output.F90.o)
..         Referenced in libgamlib.a(gioH5.F90.o)
..         Referenced in libgamlib.a(init.F90.o)'

.. A1: This may have to do with PATH, e.g., invalid directories are included in $PATH and/or $PYTHONPATH.

.. Q2. What's the code unit for current density from MHD output?

.. A2. gB0/gx0/mu0 = 4.58nT/6.38e6m/(4pi*1e-7)=5.7139e-10 A/m^2. Check the units of other variables in h5 file with h5info in Matlab or similar functions in python etc.
