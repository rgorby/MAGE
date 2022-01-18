#!/usr/bin/env python
#Generates HDF-5 helio grid
#
import kaipy.gamera.gamGrids as gg
import numpy as np

#Output name
#For inner helio
fOut = "/glade/u/home/elenap/gameraHelio/grid/heliogrid-test.h5"
#For outer helio 1-10 AU
#fOut = "/glade/u/home/elenap/gameraHelio/OHelio/grid/heliogrid-OH.h5"

#IJK = r,theta,phi

#Theta bounds (in units of 0-1)
tMin = 0.1
tMax = 0.9
#Radial bounds in R_S for inner helio
Rin = 21.5
Rout = 215.
#Radial bounds in AU for outer helio 1-10 AU
#Rin = 1.
#Rout = 10.
#Grid cells
Ni = 256
Nj = 128
Nk = 256

print("Generating gamera-helio grid ...")

X3,Y3,Z3 = gg.GenKSph(Ni=Ni,Nj=Nj,Nk=Nk,Rin=Rin,Rout=Rout,tMin=tMin,tMax=tMax)

#to generate non-uniform grid for GL cme (more fine in region 0.1-0.3 AU) 
#X3,Y3,Z3 = gg.GenKSphNonUGL(Ni=Ni,Nj=Nj,Nk=Nk,Rin=Rin,Rout=Rout,tMin=tMin,tMax=tMax)
gg.WriteGrid(X3,Y3,Z3,fOut=fOut)

