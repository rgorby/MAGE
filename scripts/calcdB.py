#!/usr/bin/env python

import sys
import argparse
from argparse import RawTextHelpFormatter
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from   astropy.time import Time
import pickle

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

# just temporary
phiStart = float(sys.argv[1])
phiEnd   = float(sys.argv[2])
altitude = float(sys.argv[3])
phiStep  = 0.5  # 0.5 degree resolution in phi

ion = remix.remix(remixFile,step)

p = np.arange(phiStart,phiEnd,phiStep)*np.pi/180. 
t = (np.arange(0,45,0.5)+0.5)*np.pi/180.    # 0.5 degree resolution in theta down to 45 colat, omitting pole
R = (6380.+altitude)/6500.           # altitude in units of Rion. Note the definitions of Re and Ri

phi,theta = np.meshgrid(p,t)
x = R*np.sin(theta)*np.cos(phi)
y = R*np.sin(theta)*np.sin(phi)
z = R*np.cos(theta)

xyz = np.array([x.ravel(),y.ravel(),z.ravel()]).T

dBr,dBtheta,dBphi = ion.dB(xyz)

print("Done computing. Saving pickles.")
pickle.dump([p,t,R,dBr,dBtheta,dBphi],open('db_%03d_%03d_%03d.pkl'%(phiStart,phiEnd,altitude),'wb'))
