#! /usr/bin/env python

import numpy as np
import os,sys,glob
from scipy import interpolate
import time
import h5py
import matplotlib.pyplot as plt

import kaipy.gamhelio.wsa2gamera.params as params
import kaipy.gamhelio.lib.wsa as wsa

import kaipy.gamera.gamGrids as gg

# Parse arguments
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('ConfigFileName',help='The name of the configuration file to use',default='startup.config')
args = parser.parse_args()

# constants
# TODO: move to kaipy/kdefs.py
# TODO: add Boltzmann constant and remove from pressure calc
mp = 1.67e-24

# Read params from config file
prm = params.params(args.ConfigFileName)
Ng=prm.NO2
gamma = prm.gamma

# remember to use the same units in gamera
B0 = prm.B0
n0 = prm.n0

# TODO move to config file
T0 = 3.44e6  #2.88e6
V0 = B0/np.sqrt(4*np.pi*mp*n0)
# note, keep temperature in K (pressure is normalized in wsa.F90)

#grid parameters
tMin = prm.tMin
tMax = prm.tMax
Rin = prm.Rin
Rout = prm.Rout
Ni = prm.Ni
Nj = prm.Nj
Nk = prm.Nk 

ffits = os.path.join(os.getenv('KAIJUHOME'),prm.wsaFile)

# Generate Helio Grid
print("Generating gamera-helio grid ...")

X3,Y3,Z3 = gg.GenKSph(Ni=Ni,Nj=Nj,Nk=Nk,Rin=Rin,Rout=Rout,tMin=tMin,tMax=tMax)
gg.WriteGrid(X3,Y3,Z3,fOut=os.path.join(prm.GridDir,prm.gameraGridFile))

print("Gamera-helio grid ready!")

# Read/normalize WSA
jd_c,phi_wsa_v,theta_wsa_v,phi_wsa_c,theta_wsa_c,bi_wsa,v_wsa,n_wsa,T_wsa = wsa.read(ffits,prm.densTempInfile,prm.normalized)

# convert the units
bi_wsa /= B0
n_wsa  /= (mp*n0)
v_wsa /= V0

#convert julian date from wsa fits into modified julian date
mjd_c = jd_c - 2400000.5


# GAMERA grid
with h5py.File(os.path.join(prm.GridDir,prm.gameraGridFile),'r') as f:
    x=f['X'][:]
    y=f['Y'][:]
    z=f['Z'][:]

xc = 0.125*(x[:-1,:-1,:-1]+x[:-1,1:,:-1]+x[:-1,:-1,1:]+x[:-1,1:,1:]
            +x[1:,:-1,:-1]+x[1:,1:,:-1]+x[1:,:-1,1:]+x[1:,1:,1:])
yc = 0.125*(y[:-1,:-1,:-1]+y[:-1,1:,:-1]+y[:-1,:-1,1:]+y[:-1,1:,1:]
            +y[1:,:-1,:-1]+y[1:,1:,:-1]+y[1:,:-1,1:]+y[1:,1:,1:])
zc = 0.125*(z[:-1,:-1,:-1]+z[:-1,1:,:-1]+z[:-1,:-1,1:]+z[:-1,1:,1:]
            +z[1:,:-1,:-1]+z[1:,1:,:-1]+z[1:,:-1,1:]+z[1:,1:,1:])

# remove the ghosts from angular dimensions
R0 = np.sqrt(x[0,0,Ng]**2+y[0,0,Ng]**2+z[0,0,Ng]**2)  # radius of the inner boundary

P = np.arctan2(y[Ng:-Ng-1,Ng:-Ng-1,:],x[Ng:-Ng-1,Ng:-Ng-1,:])
P[P<0]=P[P<0]+2*np.pi
P = P % (2*np.pi)  # sometimes the very first point may be a very
                   # small negative number, which the above call sets
                   # to 2*pi. This takes care of it.

Rc = np.sqrt(xc[Ng:-Ng,Ng:-Ng,:]**2+yc[Ng:-Ng,Ng:-Ng,:]**2+zc[Ng:-Ng,Ng:-Ng,:]**2)
Pc = np.arctan2(yc[Ng:-Ng,Ng:-Ng,:],xc[Ng:-Ng,Ng:-Ng,:])
Pc[Pc<0]=Pc[Pc<0]+2*np.pi
Tc = np.arccos(zc[Ng:-Ng,Ng:-Ng,:]/Rc)

# this is fast and better than griddata in that it nicely extrapolates boundaries:
fbi      = interpolate.RectBivariateSpline(phi_wsa_c,theta_wsa_c,bi_wsa.T,kx=1,ky=1)  
br = fbi(Pc[:,0,0],Tc[0,:,0])

# Smoothing
if not prm.gaussSmoothWidth==0:
    import astropy
    from astropy.convolution import convolve,Gaussian2DKernel

    gauss=Gaussian2DKernel(width=prm.gaussSmoothWidth)
    br   =astropy.convolution.convolve(br,gauss,boundary='extend')


# Interpolate to Gamera grid
fv      = interpolate.RectBivariateSpline(phi_wsa_c,theta_wsa_c,v_wsa.T,kx=1,ky=1)  
vr = fv(Pc[:,0,0],Tc[0,:,0])

f      = interpolate.RectBivariateSpline(phi_wsa_c,theta_wsa_c,n_wsa.T,kx=1,ky=1)  
rho = f(Pc[:,0,0],Tc[0,:,0])

# TODO: clean up, define constants above if needed
# calculate temperature from total force balance
# note, this is done AFTER interpolating br and rho to the gamera grid
temp =  1.*T0/rho + (1.**2-(br)**2)*V0**2 / 2e8/1.38 * 1.67/rho   # *****  
temp_T = temp.T

# note, redefining interpolation functions we could also
# interpolate from bi_wsa as above, but then we would have to
# smooth bk, if necessary. The way we're doing it here, bk will be
# smoothed or not, dependent on whether br has been smoothed.
# note also, this has to extrapolate
fbi = interpolate.RectBivariateSpline(Pc[:,0,0],Tc[0,:,0],br,kx=1,ky=1)
fv  = interpolate.RectBivariateSpline(Pc[:,0,0],Tc[0,:,0],vr,kx=1,ky=1)

br_kface = fbi(P[:,0,0],Tc[0,:,0])
vr_kface  = fv (P[:,0,0],Tc[0,:,0])

# Scale inside ghost region
(vr,vr_kface,rho,temp,br,br_kface) = [np.dstack(prm.NO2*[var]) for var in (vr,vr_kface,rho,temp,br,br_kface)]
rho*=(R0/Rc[0,0,:Ng])**2
br*=(R0/Rc[0,0,:Ng])**2
br_kface*=(R0/Rc[0,0,:Ng])**2

# Save to file
with h5py.File(os.path.join(prm.IbcDir,prm.gameraIbcFile),'w') as hf:
    hf.attrs["MJD"] = mjd_c
    hf.create_dataset("vr",data=vr)
    hf.create_dataset("vr_kface",data=vr_kface)
    hf.create_dataset("rho",data=rho)
    hf.create_dataset("temp",data=temp)
    hf.create_dataset("br",data=br)
    hf.create_dataset("br_kface",data=br_kface)
hf.close()
