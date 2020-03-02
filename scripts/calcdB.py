#!/usr/bin/env python

import sys
import argparse
from argparse import RawTextHelpFormatter
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from   astropy.time import Time

import kaipy.kaiH5 as kaiH5
import kaipy.remix.remix as remix


remixFile = "sigmaH.2430.mix.h5"
step = 0

#x = np.zeros((89,720,1))
#y = np.zeros((1,1,64800))
#z = x+y
#print(z.shape)
#sys.exit(0)

# Initialize the remix class
#ion = remix.remix(remixFile,args.n)
ion = remix.remix(remixFile,step)

p = np.arange(90)*np.pi/180.   # 0.5 degree resolution in phi
t = np.arange(90)*np.pi/180.    # 0.5 degree resolution in theta down to 45 colat
R = (6380.+85.)/6500.           # altitude in units of Rion. Note the definitions of Re and Ri

phi,theta = np.meshgrid(p,t)
x = R*np.sin(theta)*np.cos(phi)
y = R*np.sin(theta)*np.sin(phi)
z = R*np.cos(theta)

xyz = np.array([x.ravel(),y.ravel(),z.ravel()]).T

dBr,dBtheta,dBphi = ion.dB(xyz)

print(dBr.shape)
