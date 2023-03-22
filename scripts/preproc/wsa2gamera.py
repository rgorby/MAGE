#! /usr/bin/env python

import numpy as np
import os,sys,glob
from scipy import interpolate
import time
import h5py
import matplotlib.pyplot as plt

import kaipy.gamhelio.wsa2gamera.params as params
import kaipy.gamhelio.lib.wsa as wsa
from kaipy.kdefs import *

import kaipy.gamera.gamGrids as gg

# Parse arguments
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('ConfigFileName',help='The name of the configuration file to use',default='startup.config')
args = parser.parse_args()


# Read params from config file
prm   = params.params(args.ConfigFileName)
Ng    = prm.Nghost
gamma = prm.gamma

# Normalization parameters
# remember to use the same units in gamera
B0 = prm.B0
n0 = prm.n0
V0 = B0/np.sqrt(4*np.pi*mp*n0)
TCS = prm.TCS #Temperature in the current sheet for pressure balance calculation 
nCS = prm.nCS #Density in the current sheet for pressure balance calculation 

# Grid parameters
tMin = prm.tMin
tMax = prm.tMax
Rin  = prm.Rin
Rout = prm.Rout
Ni   = prm.Ni
Nj   = prm.Nj
Nk   = prm.Nk 

ffits = prm.wsaFile

# Generate spherical helio grid
print("Generating gamera-helio grid Ni = %d, Nj  = %d, Nk = %d " % (Ni, Nj, Nk))

X3,Y3,Z3 = gg.GenKSph(Ni=Ni,Nj=Nj,Nk=Nk,Rin=Rin,Rout=Rout,tMin=tMin,tMax=tMax)
gg.WriteGrid(X3,Y3,Z3,fOut=os.path.join(prm.GridDir,prm.gameraGridFile))

if os.path.exists(prm.gameraGridFile):
    print("Grid file heliogrid.h5 is ready!")


# Read and normalize WSA
jd_c,phi_wsa_v,theta_wsa_v,phi_wsa_c,theta_wsa_c,bi_wsa,v_wsa,n_wsa,T_wsa = wsa.read(ffits,prm.densTempInfile,prm.normalized)
bi_wsa /= B0
n_wsa  /= (mp*n0)
v_wsa /= V0
#convert julian date in the center of the WSA map into modified julian date
mjd_c = jd_c - JD2MJD


# Get GAMERA grid for further interpolation
with h5py.File(os.path.join(prm.GridDir,prm.gameraGridFile),'r') as f:
    x=f['X'][:]
    y=f['Y'][:]
    z=f['Z'][:]
# Cell centers, note order of indexes [k,j,i]
xc = 0.125*(x[:-1,:-1,:-1]+x[:-1,1:,:-1]+x[:-1,:-1,1:]+x[:-1,1:,1:]
            +x[1:,:-1,:-1]+x[1:,1:,:-1]+x[1:,:-1,1:]+x[1:,1:,1:])
yc = 0.125*(y[:-1,:-1,:-1]+y[:-1,1:,:-1]+y[:-1,:-1,1:]+y[:-1,1:,1:]
            +y[1:,:-1,:-1]+y[1:,1:,:-1]+y[1:,:-1,1:]+y[1:,1:,1:])
zc = 0.125*(z[:-1,:-1,:-1]+z[:-1,1:,:-1]+z[:-1,:-1,1:]+z[:-1,1:,1:]
            +z[1:,:-1,:-1]+z[1:,1:,:-1]+z[1:,:-1,1:]+z[1:,1:,1:])

# radius of the inner boundary
R0 = np.sqrt(x[0,0,Ng]**2+y[0,0,Ng]**2+z[0,0,Ng]**2) 

# Calculate phi and theta in physical domain (excluding ghost cells)
P = np.arctan2(y[Ng:-Ng-1,Ng:-Ng-1,:],x[Ng:-Ng-1,Ng:-Ng-1,:])
P[P<0]=P[P<0]+2*np.pi
P = P % (2*np.pi)  # sometimes the very first point may be a very
                   # small negative number, which the above call sets
                   # to 2*pi. This takes care of it.

# Calculate r, phi and theta coordinates of cell centers in physical domain (excluding ghost cells)
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

# Not interpolating temperature, but calculating from the total pressure balance
# AFTER interpolating br and rho to the gamera grid
# n_CS*k*T_CS = n*k*T + Br^2/8pi  
temp = (nCS*kbltz*TCS - (br*B0)**2/8./np.pi)/(rho*n0)/kbltz
# note, keep temperature in K (pressure is normalized in wsa.F90)

#check
#print ("Max and min of temperature in MK")
#print (np.amax(temp)*1.e-6, np.amin(temp)*1.e-6)

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
(vr,vr_kface,rho,temp,br,br_kface) = [np.dstack(Ng*[var]) for var in (vr,vr_kface,rho,temp,br,br_kface)]
rho*=(R0/Rc[0,0,:Ng])**2
br*=(R0/Rc[0,0,:Ng])**2
br_kface*=(R0/Rc[0,0,:Ng])**2

# Calculating electric field component on k_edges 
# E_theta = B_phi*Vr = - Omega*R*sin(theta)/Vr*Br * Vr = - Omega*R*sin(theta)*Br 
omega=2*np.pi/Tsolar
et_kedge = - omega*R0*np.sin(Tc[:,:,Ng-1])*br_kface[:,:,-1]

# v, rho, br are normalized, temp is in [K]
with h5py.File(os.path.join(prm.IbcDir,prm.gameraIbcFile),'w') as hf:
    hf.attrs["MJD"] = mjd_c
    hf.create_dataset("vr",data=vr)
    hf.create_dataset("vr_kface",data=vr_kface)
    hf.create_dataset("rho",data=rho)
    hf.create_dataset("temp",data=temp)
    hf.create_dataset("br",data=br)
    hf.create_dataset("br_kface",data=br_kface)
    #hf.create_dataset("et_kedge",data=et_kedge)
hf.close()
