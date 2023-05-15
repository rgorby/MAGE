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

TCS = prm.TCS # Temperature in the current sheet for pressure balance calculation 
nCS = prm.nCS # Density in the current sheet for pressure balance calculation 

# Grid parameters
tMin = prm.tMin
tMax = prm.tMax
Rin  = prm.Rin
Rout = prm.Rout
Ni   = prm.Ni
Nj   = prm.Nj
Nk   = prm.Nk 

#conversions from wsa to gamera units
cms2kms = 1.e-5 # cm/s => km/s
Gs2nT = 1.e5    # Gs => nT
# Conversion for E field  1 statV/cm = 3.e7 mV/m
eScl = 3.e7

ffits = prm.wsaFile

# Generate spherical helio grid
print("Generating gamera-helio grid Ni = %d, Nj  = %d, Nk = %d " % (Ni, Nj, Nk))

X3,Y3,Z3 = gg.GenKSph(Ni=Ni,Nj=Nj,Nk=Nk,Rin=Rin,Rout=Rout,tMin=tMin,tMax=tMax)
gg.WriteGrid(X3,Y3,Z3,fOut=os.path.join(prm.GridDir,prm.gameraGridFile))

if os.path.exists(prm.gameraGridFile):
    print("Grid file heliogrid.h5 is ready!")

# Read WSA
jd_c,phi_wsa_v,theta_wsa_v,phi_wsa_c,theta_wsa_c,bi_wsa,v_wsa,n_wsa,T_wsa = wsa.read(ffits,prm.densTempInfile,prm.normalized)
# Units of WSA input
# bi_wsa in [Gs] 
# v_wsa in [cm/s]
# n_wsa in [g cm-3]
# T_wsa in [K]

# convert julian date in the center of the WSA map into modified julian date
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

r = np.sqrt(x**2+y**2+z**2)

# Calculate phi and theta in physical domain (excluding ghost cells)
P = np.arctan2(y[Ng:-Ng,Ng:-Ng,:],x[Ng:-Ng,Ng:-Ng,:])
P[P<0]=P[P<0]+2*np.pi
#P = P % (2*np.pi)  # sometimes the very first point may be a very
                   # small negative number, which the above call sets
                   # to 2*pi. This takes care of it.
T = np.arccos(z[Ng:-Ng,Ng:-Ng,:]/r[Ng:-Ng,Ng:-Ng,:])

#grid for inner i-ghost region; output to innerbc.h5
P_out = P[:,:,0:Ng+1]
T_out = T[:,:,0:Ng+1]
R_out = r[Ng:-Ng,Ng:-Ng,0:Ng+1]

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
temp = (nCS*kbltz*TCS - br**2/8./np.pi)*Mp_cgs/rho/kbltz
#temperature in [K]

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

br_kface = fbi(P[:-1,0,0],Tc[0,:,0]) #(Nk,Nj)
vr_kface  = fv (P[:-1,0,0],Tc[0,:,0]) #(Nk,Nj)

# before applying scaling inside ghost region
# get br values to the left of an edge for E_theta calculation
br_kedge = np.roll(br,1, axis=1)

# Scale inside ghost region
(vr,vr_kface,rho,temp,br,br_kface) = [np.dstack(Ng*[var]) for var in (vr,vr_kface,rho,temp,br,br_kface)]
rho*=(R0/Rc[0,0,:Ng])**2
br*=(R0/Rc[0,0,:Ng])**2
br_kface*=(R0/Rc[0,0,:Ng])**2

# Calculating E-field component on k_edges in [mV/m]
# E_theta = B_phi*Vr/c = - Omega*R*sin(theta)/Vr*Br * Vr/c = - Omega*R*sin(theta)*Br/c
omega = 2*np.pi/(Tsolar*Day2s) # [1/s]
# Theta at centers of k-faces (== theta at kedges)
Tcf = 0.25*(T[:,:-1,:-1] + T[:,1:,1:] + T[:,:-1,1:] + T[:,1:,:-1])
et_kedge = - omega*R0*Rsolar*np.sin(Tcf[:-1,:,Ng-1])*br_kedge/vc_cgs #[statV/cm]

# Unit conversion agreement. Input to GAMERA innerbc.h5 has units V[km/s], Rho[cm-3], T[K], B[nT], E[mV/m]
vr *= cms2kms
vr_kface *= cms2kms
rho /= Mp_cgs
br *= Gs2nT
br_kface *= Gs2nT
et_kedge *= eScl

with h5py.File(os.path.join(prm.IbcDir,prm.gameraIbcFile),'w') as hf:
    hf.create_dataset("X", data=P_out)
    hf.create_dataset("Y", data=T_out)
    hf.create_dataset("Z", data=R_out)
    grname = "Step#0"
    grp = hf.create_group(grname)
    grp.attrs.create("MJD", mjd_c)
    grp.create_dataset("vr",data=vr)
    grp.create_dataset("vr_kface",data=vr_kface) # size (Nk,Nj,Ng)
    grp.create_dataset("rho",data=rho)
    grp.create_dataset("temp",data=temp)
    grp.create_dataset("br",data=br)
    grp.create_dataset("br_kface",data=br_kface) # size (Nk,Nj,Ng)
    grp.create_dataset("et_kedge",data=et_kedge)  # size (Nk, Nj)
hf.close()
