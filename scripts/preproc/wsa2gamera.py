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

#----------- PARSE ARGUMENTS ---------#
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('ConfigFileName',help='The name of the configuration file to use',default='startup.config')
args = parser.parse_args()
#----------- PARSE ARGUMENTS ---------#

# Read params from config file
prm = params.params(args.ConfigFileName)
Ng=prm.NO2
gamma = prm.gamma
B0 = prm.B0
n0 = prm.n0
T0 = 3.44e6  #2.88e6

#grid parameters
tMin = prm.tMin
tMax = prm.tMax
Rin = prm.Rin
Rout = prm.Rout
Ni = prm.Ni
Nj = prm.Nj
Nk = prm.Nk 

# constants
mp = 1.67e-24

#----------GENERATE HELIO GRID------

print("Generating gamera-helio grid ...")

X3,Y3,Z3 = gg.GenKSph(Ni=Ni,Nj=Nj,Nk=Nk,Rin=Rin,Rout=Rout,tMin=tMin,tMax=tMax)

#to generate non-uniform grid for GL cme (more fine in region 0.1-0.3 AU) 
#X3,Y3,Z3 = gg.GenKSphNonUGL(Ni=Ni,Nj=Nj,Nk=Nk,Rin=Rin,Rout=Rout,tMin=tMin,tMax=tMax)
gg.WriteGrid(X3,Y3,Z3,fOut=os.path.join(prm.GridDir,prm.gameraGridFile))

print("Gamera-helio grid ready!")

#----------GENERATE HELIO GRID------


############### WSA STUFF #####################
jd_c,phi_wsa_v,theta_wsa_v,phi_wsa_c,theta_wsa_c,bi_wsa,v_wsa,n_wsa,T_wsa = wsa.read(prm.wsaFile,prm.densTempInfile,prm.normalized)

# convert the units; remember to use the same units in gamera
# TODO: probably store units in the h5 file?
# B0   = 1.e-3 Gs
# n0   = 200./cc

V0 = B0/np.sqrt(4*np.pi*mp*n0)

bi_wsa /= B0
n_wsa  /= (mp*n0)
v_wsa /= V0
#convert julian date from wsa fits into modified julian date
mjd_c = jd_c - 2400000.5
# keep temperature in K
############### WSA STUFF #####################


############### GAMERA STUFF #####################

# GAMERA GRID
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

############### SMOOTHING #####################
if not prm.gaussSmoothWidth==0:
    import astropy
    from astropy.convolution import convolve,Gaussian2DKernel

    gauss=Gaussian2DKernel(width=prm.gaussSmoothWidth)
    br   =astropy.convolution.convolve(br,gauss,boundary='extend')


############### INTERPOLATE AND DUMP #####################
fv      = interpolate.RectBivariateSpline(phi_wsa_c,theta_wsa_c,v_wsa.T,kx=1,ky=1)  
vr = fv(Pc[:,0,0],Tc[0,:,0])

f      = interpolate.RectBivariateSpline(phi_wsa_c,theta_wsa_c,n_wsa.T,kx=1,ky=1)  
rho = f(Pc[:,0,0],Tc[0,:,0])

#f      = interpolate.RectBivariateSpline(phi_wsa_c,theta_wsa_c,T_wsa.T,kx=1,ky=1)  
#temp = f(Pc[:,0,0],Tc[0,:,0])
temp =  1.*T0/rho + (1.**2-(br)**2)*V0**2 / 2e8/1.38 * 1.67/rho   # *****  
temp_T = temp.T

pressure =   ((br)**2)*V0**2 /2.*mp*n0 *0.1   + (n0*rho* temp)*1.38e-16 *0.1
pressure_therm = (n0*rho* temp)*1.38e-16 * 0.1
pressure_B = (br)**2 *V0**2 / 2.*mp*n0 *0.1


# Plot different variables 
#nk = 256  # 1024
#nj = 128   # 512

#Pcc   = np.arange(nk)*360./nk
#Tcc = np.arange(nj)*0.8*180./nj + 0.1*180.

#fig = plt.figure(figsize=(7.5,3))
#plt.pcolormesh(Pcc,Tcc, temp[:,::-1].T,cmap=cm.RdBu,rasterized=True)
#plt.xlabel(r'$\Phi$ [$^\circ$]')
#plt.ylabel(r'$\Theta$ [$^\circ$]')
#plt.colorbar()
#plt.axes().set_aspect('equal')
#plt.savefig('./temp_wsa2gam.pdf',dpi=300)

#fig = plt.figure(figsize=(7.5,3))
#plt.pcolormesh(Pcc,Tcc, pressure_therm[:,::-1].T,cmap=cm.plasma,vmin=pressure.min()*0.99,vmax=pressure.max()*1.01,rasterized=True)
#plt.xlabel(r'$\Phi$ [$^\circ$]')
#plt.ylabel(r'$\Theta$ [$^\circ$]')
#plt.colorbar()
#plt.axes().set_aspect('equal')
#plt.savefig('./thermal_Pres_wsa2gam.pdf',dpi=300)

#fig = plt.figure(figsize=(7.5,3))
#plt.pcolormesh(Pcc,Tcc, pressure[:,::-1].T,cmap=cm.plasma,vmin=pressure.min()*0.9,vmax=pressure.max()*1.1,rasterized=True)
#plt.xlabel(r'$\Phi$ [$^\circ$]')
#plt.ylabel(r'$\Theta$ [$^\circ$]')
#plt.colorbar()
#plt.axes().set_aspect('equal')
#plt.savefig('./total_P_wsa2gam.pdf',dpi=300)
#fig = plt.figure(figsize=(7.5,3))
#plt.pcolormesh(Pcc,Tcc, pressure_B[:,::-1].T,cmap=cm.plasma,vmin=pressure_B.min()*0.9,vmax=pressure_B.max()*1.1,rasterized=True)
#plt.xlabel(r'$\Phi$ [$^\circ$]')
#plt.ylabel(r'$\Theta$ [$^\circ$]')
#plt.colorbar()
#plt.axes().set_aspect('equal')
#plt.savefig('./magnetic_P_wsa2gam.pdf',dpi=300)

#fig = plt.figure(figsize=(7.5,3))
#plt.pcolormesh(Pcc,Tcc, br[:,::-1].T,cmap=cm.RdBu,vmin=br.min(),vmax=br.max(),rasterized=True)
#plt.xlabel(r'$\Phi$ [$^\circ$]')
#plt.ylabel(r'$\Theta$ [$^\circ$]')
#plt.colorbar()
#plt.axes().set_aspect('equal')
#plt.savefig('./Br_wsa2gam.pdf',dpi=300)

#fig = plt.figure(figsize=(7.5,3))
#plt.pcolormesh(Pcc,Tcc, vr[:,::-1].T,cmap=cm.plasma,vmin=vr.min(),vmax=vr.max(),rasterized=True)
#plt.xlabel(r'$\Phi$ [$^\circ$]')
#plt.ylabel(r'$\Theta$ [$^\circ$]')
#plt.colorbar()
#plt.axes().set_aspect('equal')
#plt.savefig('./Vr_wsa2gam.pdf',dpi=300)

#fig = plt.figure(figsize=(7.5,3))
#plt.pcolormesh(Pcc,Tcc, rho[:,::-1].T,cmap=cm.RdBu,vmin=rho.min(),vmax=rho.max(),rasterized=True)
#plt.xlabel(r'$\Phi$ [$^\circ$]')
#plt.ylabel(r'$\Theta$ [$^\circ$]')
#plt.colorbar()
#plt.axes().set_aspect('equal')
#plt.savefig('./Density_wsa2gam.pdf',dpi=300)

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

with h5py.File(os.path.join(prm.IbcDir,prm.gameraIbcFile),'w') as hf:
    hf.attrs["MJD"] = mjd_c
    hf.create_dataset("vr",data=vr)
    hf.create_dataset("vr_kface",data=vr_kface)
    hf.create_dataset("rho",data=rho)
    hf.create_dataset("temp",data=temp)
    hf.create_dataset("br",data=br)
    hf.create_dataset("br_kface",data=br_kface)
hf.close()
