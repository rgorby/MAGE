#!/usr/bin/env python

import h5py
import numpy as np
import sys

TW = 1.4493e+5     #Default temperature, K => 0.01 nPa
nW = 5          #Default density, #/cc
VxW = 400.0     #Default wind, km/s
f107val = 100.0 #Default f10.7 flux
tilt = 0.0      #Default dipole tilt, radians
mjd0 = 58767.0  #Default MJD, set for 2019-10-11 00:00:00                                                                                                                           
Bx0 = 0.0 #Default Bx offset for planar front, keep at zero
ByC = 0.0 #Default By coefficient used to calculate Bx, include if want tilted field
BzC = 0.0 #Default Bz coefficient used to calculate Bx, include if want tilted field                                                                                                                
fOut = "bcwind_2.h5"

#Time bounds [hours]                                                                                                                                                                                 
tMin = 0.0
tMax = 16.0
dt = 15.0 #Cadence [s]                                                                                                                                                                               

SimT = (tMax-tMin)*60.0*60.0
NumT = np.int( np.ceil(SimT/dt)+1 )

print("Generating %d slices, T=[%5.2f,%5.2f]"%(NumT,tMin,tMax))

T = np.linspace(tMin,tMax,NumT)
D = np.zeros(NumT)
Temp = np.zeros(NumT)
Vx = np.zeros(NumT)
Vy = np.zeros(NumT)
Vz = np.zeros(NumT)
Bx = np.zeros(NumT)
By = np.zeros(NumT)
Bz = np.zeros(NumT)
f107 = np.zeros(NumT)
ThT = np.zeros(NumT)
mjd = np.zeros(NumT)
symh = np.zeros(NumT)

tWin = 1.0 #Window times [hr]                                                                                                                                                                        
for i in range(NumT):
    t = T[i] #Time in hours                                                                                                                                                                          
    if (t <= tWin):
        D[i] = nW
        Vx[i] = -VxW
        Temp[i] = TW
        f107[i] = f107val
        ThT[i] = tilt
        mjd[i] = mjd0 + T[i]/24.0
    elif (t <= 3.5*tWin):
        D[i] = nW
        Vx[i] = -VxW
        Temp[i] = TW
        Bz[i] = 5.0
        f107[i] = f107val
        ThT[i] = tilt
        mjd[i] = mjd0 + T[i]/24.0
    elif (t <= 6.0*tWin):
        D[i] = nW
        Vx[i] = -VxW
        Temp[i] = TW
        Bz[i] = -5.0
        f107[i] = f107val
        ThT[i] = tilt
        mjd[i] = mjd0 + T[i]/24.0
    elif (t <= 8.0*tWin):
        D[i] = nW
        Vx[i] = -VxW
        Temp[i] = TW
        Bz[i] = 5.0
        f107[i] = f107val
        ThT[i] = tilt
        mjd[i] = mjd0 + T[i]/24.0
    else:
        D[i] = 2.0*nW
        Vx[i] = -VxW
        Temp[i] = TW
        Bz[i] = -15.0
        f107[i] = f107val
        ThT[i] = tilt
        mjd[i] = mjd0 + T[i]/24.0

#Write solar wind                                                                                                                                                                                    
#t,D,V,Temp,B = [s],[#/cm3],[m/s],[K],[nT]                                                                                                                                                            

oTScl = (60*60.0) #hr->s                                                                                                                                                                             
oDScl = 1.0
oVScl = 1.0e+3 #km/s->m/s                                                                                                                                                                            
oTempScl = 1.0
oBScl = 1.0


with h5py.File(fOut,'w') as hf:
    hf.create_dataset("T" ,data=oTScl*T)


    hf.create_dataset("symh" ,data=symh)
    hf.create_dataset("D" ,data=oDScl*D)
    hf.create_dataset("Temp" ,data=oTempScl*Temp)
    hf.create_dataset("Vx",data=oVScl*Vx)
    hf.create_dataset("Vy",data=oVScl*Vy)
    hf.create_dataset("Vz",data=oVScl*Vz)
    hf.create_dataset("Bx",data=oBScl*Bx)
    hf.create_dataset("By",data=oBScl*By)
    hf.create_dataset("Bz",data=oBScl*Bz)
    hf.create_dataset("tilt",data=ThT)
    hf.create_dataset("f10.7",data=f107)
    hf.create_dataset("MJD",data=mjd)
    hf.create_dataset("Bx0",data=Bx0)
    hf.create_dataset("ByC",data=ByC)
    hf.create_dataset("BzC",data=BzC)

    
