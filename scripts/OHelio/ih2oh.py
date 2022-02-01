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
# constants
mp = 1.67e-24
kb = 1.38e-16

#grid parameters
tMin = prm.tMin
tMax = prm.tMax
Rin = prm.Rin
Rout = prm.Rout
Ni = prm.Ni
Nj = prm.Nj
Nk = prm.Nk

#normalization in IH
B0 = prm.B0
n0 = prm.n0
V0 = B0/np.sqrt(4*np.pi*mp*n0)
T0 = B0*B0/4/np.pi/n0/kb #in K p = nkT

print ("inner helio normalization")
print (B0, n0, V0, T0)

#normalization in OH
B0OH = 5.e-5 # [Gs] 5 nT = 5.e-5 Gs
n0OH = 10 # [cm-3]
V0OH = B0OH/np.sqrt(4*np.pi*mp*n0OH) #Alfven speed at 1 AU 34.5 [km/s]
T0OH = B0OH*B0OH/4/np.pi/n0OH/kb #in K p = nkT

print ("outer helio units")
print (B0OH, n0OH, V0OH, T0OH)

#----------GENERATE HELIO GRID------

print("Generating gamera-Ohelio grid ...")

X3,Y3,Z3 = gg.GenKSph(Ni=Ni,Nj=Nj,Nk=Nk,Rin=Rin,Rout=Rout,tMin=tMin,tMax=tMax)

#to generate non-uniform grid for GL cme (more fine in region 0.1-0.3 AU) 
#X3,Y3,Z3 = gg.GenKSphNonUGL(Ni=Ni,Nj=Nj,Nk=Nk,Rin=Rin,Rout=Rout,tMin=tMin,tMax=tMax)
gg.WriteGrid(X3,Y3,Z3,fOut=os.path.join(prm.GridDir,prm.gameraGridFile))

print("Gamera-Ohelio grid ready!")

#----------GENERATE HELIO GRID------


############### READ GAMERA solution at 1 AU #####################
f = h5py.File(prm.wsaFile,'r')
#the latest Step saved in inner helio solution wsa.h5
step = 'Step#2'

f[step].attrs.keys()
Nphi, Nth, Nr = np.shape(f[step]['Vx'])

#coordinates of cell centers
x = 0.125*(f['X'][:-1,:-1,:-1]+f['X'][:-1,:-1,1:]+f['X'][:-1,1:,:-1]+f['X'][:-1,1:,1:]+
        f['X'][1:,:-1,:-1]+f['X'][1:,:-1,1:]+f['X'][1:,1:,:-1]+f['X'][1:,1:,1:])
y = 0.125*(f['Y'][:-1,:-1,:-1]+f['Y'][:-1,:-1,1:]+f['Y'][:-1,1:,:-1]+f['Y'][:-1,1:,1:]+
        f['Y'][1:,:-1,:-1]+f['Y'][1:,:-1,1:]+f['Y'][1:,1:,:-1]+f['Y'][1:,1:,1:])
z = 0.125*(f['Z'][:-1,:-1,:-1]+f['Z'][:-1,:-1,1:]+f['Z'][:-1,1:,:-1]+f['Z'][:-1,1:,1:]+
        f['Z'][1:,:-1,:-1]+f['Z'][1:,:-1,1:]+f['Z'][1:,1:,:-1]+f['Z'][1:,1:,1:])

r = np.sqrt(x[:]**2 + y[:]**2 + z[:]**2)
rxy = np.sqrt(x[:]**2 + y[:]**2)

#phi and theta of centers
theta = np.arccos(z/r)
phi = np.arctan2(y[:], x[:])
phi[phi<0]=phi[phi<0]+2*np.pi

theta_wsa_c = theta[0,:,0]
phi_wsa_c = phi[:,0,0]

print ("grid dimensions from 1 AU input solution")
print (theta_wsa_c.shape, phi_wsa_c.shape)

#these are normilized according to inner helio normalization
Vr = (f[step]['Vx'][:]*x[:] + f[step]['Vy'][:]*y[:] + f[step]['Vz'][:]*z[:])/r[:]
#Br = f[step]['Br'][:]
Br = (f[step]['Bx'][:]*x[:] + f[step]['By'][:]*y[:] + f[step]['Bz'][:]*z[:])/r[:]
Rho = f[step]['D'][:]
T = f[step]['P'][:]/f[step]['D'][:]

#take solution from the last cell in i, already normilized
#use wsa variable names for now
bi_wsa = Br[:,:,Nr-1]
v_wsa = Vr[:,:,Nr-1]
n_wsa = Rho[:,:,Nr-1]
T_wsa = T[:,:,Nr-1]

print ("1AU arrays")
print (bi_wsa.shape, v_wsa.shape, n_wsa.shape, T_wsa.shape)

#renormalize inner helio solution
bi_wsa = bi_wsa * B0/B0OH
n_wsa = n_wsa * n0/n0OH
v_wsa = v_wsa * V0/V0OH
T_wsa = T_wsa * T0     
# keep temperature in K

#######Interpolate to GAMERA grid###########
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

#cell corners including ghost cells
r = np.sqrt(x[:]**2+y[:]**2+z[:]**2)

#corners of physical cells
P = np.arctan2(y[Ng:-Ng,Ng:-Ng,:],x[Ng:-Ng,Ng:-Ng,:])
P[P<0] += 2*np.pi
#P = P % (2*np.pi)  # sometimes the very first point may be a very
                   # small negative number, which the above call sets
                   # to 2*pi. This takes care of it.
T = np.arccos(z[Ng:-Ng,Ng:-Ng,:]/r[Ng:-Ng,Ng:-Ng,:])

#grid (corners) for output into innerbc.h5
P_out = P[:,:,0:Ng+1]
T_out = T[:,:,0:Ng+1]
R_out = r[Ng:-Ng,Ng:-Ng,0:Ng+1]
print ("shapes of output phi and theta ", P_out.shape, T_out.shape, R_out.shape)

#centers
Rc = np.sqrt(xc[Ng:-Ng,Ng:-Ng,:]**2+yc[Ng:-Ng,Ng:-Ng,:]**2+zc[Ng:-Ng,Ng:-Ng,:]**2)
Pc = np.arctan2(yc[Ng:-Ng,Ng:-Ng,:],xc[Ng:-Ng,Ng:-Ng,:])
Pc[Pc<0] += 2*np.pi
Tc = np.arccos(zc[Ng:-Ng,Ng:-Ng,:]/Rc)


# this is fast and better than griddata in that it nicely extrapolates boundaries:
fbi      = interpolate.RectBivariateSpline(phi_wsa_c,theta_wsa_c,bi_wsa,kx=1,ky=1)  
br = fbi(Pc[:,0,0],Tc[0,:,0])

############### SMOOTHING #####################
if not prm.gaussSmoothWidth==0:
    import astropy
    from astropy.convolution import convolve,Gaussian2DKernel

    gauss=Gaussian2DKernel(width=prm.gaussSmoothWidth)
    br   =astropy.convolution.convolve(br,gauss,boundary='extend')


############### INTERPOLATE AND DUMP #####################
fv      = interpolate.RectBivariateSpline(phi_wsa_c,theta_wsa_c,v_wsa,kx=1,ky=1)  
vr = fv(Pc[:,0,0],Tc[0,:,0])

f      = interpolate.RectBivariateSpline(phi_wsa_c,theta_wsa_c,n_wsa,kx=1,ky=1)  
rho = f(Pc[:,0,0],Tc[0,:,0])

fT      = interpolate.RectBivariateSpline(phi_wsa_c,theta_wsa_c,T_wsa,kx=1,ky=1)  
temp = fT(Pc[:,0,0],Tc[0,:,0])

fbi = interpolate.RectBivariateSpline(Pc[:,0,0],Tc[0,:,0],br,kx=1,ky=1)
fv  = interpolate.RectBivariateSpline(Pc[:,0,0],Tc[0,:,0],vr,kx=1,ky=1)

br_kface = fbi(P[:,0,0],Tc[0,:,0])
vr_kface  = fv (P[:,0,0],Tc[0,:,0])

# Scale inside ghost region
(vr,vr_kface,rho,temp,br,br_kface) = [np.dstack(prm.NO2*[var]) for var in (vr,vr_kface,rho,temp,br,br_kface)]
rho *= (R0/Rc[0,0,:Ng])**2
br *= (R0/Rc[0,0,:Ng])**2
br_kface *= (R0/Rc[0,0,:Ng])**2

#FIX: 
#For now I use wsa.h5 which do not have mjd inside
#so hardcoded using WSA fits value + 200*4637/60/60/24 [number of days]
mjd_c = 58005.83415

print ("writing out innerbc.h5...")
with h5py.File(os.path.join(prm.IbcDir,prm.gameraIbcFile),'w') as hf:
    hf.attrs["MJD"] = mjd_c
    hf.create_dataset("vr",data=vr)
    hf.create_dataset("vr_kface",data=vr_kface)
    hf.create_dataset("rho",data=rho)
    hf.create_dataset("temp",data=temp)
    hf.create_dataset("br",data=br)
    hf.create_dataset("br_kface",data=br_kface)
hf.close()

#innerbc to plot in Paraview
with h5py.File(os.path.join(prm.IbcDir,'innerbc_OHighostgr.h5'),'w') as hfg:
    hfg.create_dataset("X", data=P_out)
    hfg.create_dataset("Y", data=T_out)
    hfg.create_dataset("Z", data=R_out)
    grname = "Step#0"
    grp = hfg.create_group(grname)
    grp.attrs.create("MJD", mjd_c)
    #grp.attrs.create("time", time_sec)
    grp.create_dataset("vr",data=vr)
    grp.create_dataset("vr_kface",data=vr_kface)
    grp.create_dataset("rho",data=rho)
    grp.create_dataset("temp",data=temp)
    grp.create_dataset("br",data=br)
    grp.create_dataset("br_kface",data=br_kface)
hfg.close
