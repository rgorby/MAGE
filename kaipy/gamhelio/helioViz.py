#Various helper routines for heliosphere quick look plots

import argparse
from argparse import RawTextHelpFormatter
import matplotlib as mpl
import numpy as np

import kaipy.kaiViz as kv
import kaipy.gamhelio.heliosphere as hsph
from kaipy.kdefs import *
import os


VMax = 800.
VMin = 300.
MagVCM = "inferno"
#MagVCM = "rainbow"

#inner helio
DMax = 150.
DMin = 2000.
DCM = "copper_r"

#limits for iSlice
#21.5 R_S
#D0Max = 1000.
#D0Min = 300.
#1 au
D0Max = 1.
D0Min = 15.
D0CM = "copper_r"

TMin = 0.2
TMax = 2.
TCM = "copper"

T0Min = 0.01
T0Max = 0.25

BMax = 150.
BMin = -150.
#BMax = 5.
#BMin = -5.
BCM = "coolwarm"

B0Min = -4.
B0Max = 4.

colorProf = "tab:orange"
#Function to Add different size options to argument
#not used for helio right now
def AddSizeArgs(parser):
	parser.add_argument('-small' , action='store_true', default=False,help="Use smaller domain bounds (default: %(default)s)")
	parser.add_argument('-big'   , action='store_true', default=False,help="Use larger domain bounds (default: %(default)s)")
	parser.add_argument('-bigger', action='store_true', default=False,help="Use larger-er domain bounds (default: %(default)s)")
	parser.add_argument('-huge'  , action='store_true', default=False,help="Use huge domain bounds (default: %(default)s)")

#Return domain size from parsed arguments; see msphViz for options
def GetSizeBds(pic):
	if (pic == "pic1" or pic == "pic2"):
                #for inner helio
		xyBds = [-220.,220.,-220.,220.]
                #for 1-10 au helio
                #xyBds = [-10.,10.,-10.,10.]
	elif (pic == "pic3"):
		xyBds = [0.,360.,-75.,75.]
	elif (pic == "pic4"):
                xyBds = [0.,360.,-90.,90.]
	elif (pic == "pic5"):
		xyBds = [20.,220.,1.,2000.]
	else:		
		print ("No pic type specified.")
	return xyBds

#Plot speed in equatorial plane
def PlotEqMagV(gsph,nStp,xyBds,Ax,AxCB=None,doClear=True,doDeco=True):
	vMagV = kv.genNorm(VMin, VMax, doLog=False, midP=None)
	
	if (AxCB is not None):
		#Add the colorbar to AxCB
		AxCB.clear()
		kv.genCB(AxCB,vMagV,"Speed [km/s]",cM=MagVCM,Ntk=7)

	#Now do main plotting
	if (doClear):
		Ax.clear()

	MagV = gsph.eqMagV(nStp)
	Ax.pcolormesh(gsph.xxi,gsph.yyi,MagV,cmap=MagVCM,norm=vMagV)

	kv.SetAx(xyBds,Ax)

	if (doDeco):
		Ax.set_xlabel('X [R_S]')
		Ax.set_ylabel('Y [R_S]')
	return MagV

#Plot speed in meridional plane Y=0
def PlotMerMagV(gsph,nStp,xyBds,Ax,AxCB=None,doClear=True,doDeco=True):
        vMagV = kv.genNorm(VMin, VMax, doLog=False, midP=None)

        if (AxCB is not None):
                #Add the colorbar to AxCB
                AxCB.clear()
                kv.genCB(AxCB,vMagV,"Speed [km/s]",cM=MagVCM,Ntk=7)

        if (doClear):
                Ax.clear()
	#r is for +X plane and l is for -X plane
        Vr, Vl = gsph.MerMagV(nStp)
	#cell corners
        xr, zr, xl, zl, r = gsph.MeridGridHalfs()
        Ax.pcolormesh(xr,zr,Vr,cmap=MagVCM,norm=vMagV)
        Ax.pcolormesh(xl,zl,Vl,cmap=MagVCM,norm=vMagV)

        kv.SetAx(xyBds,Ax)

        if (doDeco):
                Ax.set_xlabel('X [R_S]')
                Ax.set_ylabel('Z [R_S]')
        return Vr, Vl

#Plot normalized density n(r/r0)^2 in meridional plane Y=0
def PlotMerDNorm(gsph,nStp,xyBds,Ax,AxCB=None,doClear=True,doDeco=True):
	vD = kv.genNorm(DMin, DMax, doLog=False, midP=None)

	if (AxCB is not None):
		#Add the colorbar to AxCB
		AxCB.clear()
		kv.genCB(AxCB,vD,r"Density n$(r/r_0)^2$ [cm$^{-3}$]",cM=DCM,Ntk=7)

	if (doClear):
		Ax.clear()

	Dr, Dl = gsph.MerDNrm(nStp)
	xr, zr, xl, zl, r = gsph.MeridGridHalfs()
	Ax.pcolormesh(xr,zr,Dr,cmap=DCM,norm=vD, shading='auto')
	Ax.pcolormesh(xl,zl,Dl,cmap=DCM,norm=vD, shading='auto')

	kv.SetAx(xyBds,Ax)

	if (doDeco):
		Ax.set_xlabel('X [R_S]')
		Ax.set_ylabel('Z [R_S]')
		Ax.yaxis.tick_right()
		Ax.yaxis.set_label_position('right')
	return Dr, Dl

#Plot normalized Br Br(r/r0)^2 in meridional plane Y=0
def PlotMerBrNorm(gsph,nStp,xyBds,Ax,AxCB=None,doClear=True,doDeco=True):
	vB = kv.genNorm(BMin, BMax, doLog=False, midP=None)

	if (AxCB is not None):
		#Add the colorbar to AxCB
		AxCB.clear()
		kv.genCB(AxCB,vB,r'Radial MF B$_r$(r/r$_0)^2$ [nT]',cM=BCM,Ntk=7)
	if (doClear):
		Ax.clear()

	Br_r, Br_l = gsph.MerBrNrm(nStp)
	xr, zr, xl, zl, r = gsph.MeridGridHalfs()
	Ax.pcolormesh(xr,zr,Br_r,cmap=BCM,norm=vB,shading='auto')
	Ax.pcolormesh(xl,zl,Br_l,cmap=BCM,norm=vB,shading='auto')
	#plot heliospheric current sheet
	#cell-cent coords first
	xr_c = 0.25*( xr[:-1,:-1]+xr[:-1,1:]+xr[1:,:-1]+xr[1:,1:] )
	zr_c = 0.25*( zr[:-1,:-1]+zr[:-1,1:]+zr[1:,:-1]+zr[1:,1:] )
	xl_c = 0.25*( xl[:-1,:-1]+xl[:-1,1:]+xl[1:,:-1]+xl[1:,1:] )
	zl_c = 0.25*( zl[:-1,:-1]+zl[:-1,1:]+zl[1:,:-1]+zl[1:,1:] )
	#plot Br=0
	Ax.contour(xr_c,zr_c,Br_r,[0.],colors='black')
	Ax.contour(xl_c,zl_c,Br_l,[0.],colors='black')
	kv.SetAx(xyBds,Ax)

	if (doDeco):
		Ax.set_xlabel('X [R_S]')
		Ax.set_ylabel('Z [R_S]')
		Ax.yaxis.tick_right()
		Ax.yaxis.set_label_position('right')
	return Br_r, Br_l

#Plot normalized temperature T(r/r0) in meridional plane
def PlotMerTemp(gsph,nStp,xyBds,Ax,AxCB=None,doClear=True,doDeco=True):
	vT = kv.genNorm(TMin, TMax, doLog=False, midP=None)

	if (AxCB is not None):
		AxCB.clear()
		kv.genCB(AxCB,vT, r'Temperature T(r/r$_0$) [MK]',cM=TCM,Ntk=7)
	if (doClear):
		Ax.clear()

	Tempr, Templ = gsph.MerTemp(nStp)
	xr, zr, xl, zl, r = gsph.MeridGridHalfs()
	Ax.pcolormesh(xr,zr,Tempr,cmap=TCM,norm=vT)
	Ax.pcolormesh(xl,zl,Templ,cmap=TCM,norm=vT)

	kv.SetAx(xyBds,Ax)

	if (doDeco):
		Ax.set_xlabel('X [R_S]')
		Ax.set_ylabel('Z [R_S]')
	return Tempr, Templ

#Plot normalized density in equatorial plane n(r/r0)^2
def PlotEqD(gsph,nStp,xyBds,Ax,AxCB=None,doClear=True,doDeco=True):
	vD = kv.genNorm(DMin, DMax, doLog=False, midP=None)
	
	if (AxCB is not None):
		#Add the colorbar to AxCB
		AxCB.clear()
		kv.genCB(AxCB,vD,r"Density n(r/r$_0)^2$ [cm$^{-3}$]",cM=DCM,Ntk=7)

	#Now do main plotting
	if (doClear):
		Ax.clear()

	NormD = gsph.eqNormD(nStp)
	Ax.pcolormesh(gsph.xxi,gsph.yyi,NormD,cmap=DCM,norm=vD)

	kv.SetAx(xyBds,Ax)

	if (doDeco):
		Ax.set_xlabel('X [R_S]')
		Ax.set_ylabel('Y [R_S]')
		Ax.yaxis.tick_right()
		Ax.yaxis.set_label_position('right')
	return NormD

#Plot normalized Temperature in equatorial plane T(r/r0)
def PlotEqTemp(gsph,nStp,xyBds,Ax,AxCB=None,doClear=True,doDeco=True):
	vT = kv.genNorm(TMin, TMax, doLog=False, midP=None)

	if (AxCB is not None):
		AxCB.clear()
		kv.genCB(AxCB,vT,r"Temperature T(r/r$_0$) [MK]",cM=TCM,Ntk=7)
	if (doClear):
		Ax.clear()

	Temp = gsph.eqTemp(nStp)
	Ax.pcolormesh(gsph.xxi,gsph.yyi,Temp,cmap=TCM,norm=vT)
	
	kv.SetAx(xyBds,Ax)

	if (doDeco):
		Ax.set_xlabel('X [R_S]')
		Ax.set_ylabel('Y [R_S]')
	return Temp

#Plor Br in equatorial plane
def PlotEqBr(gsph,nStp,xyBds,Ax,AxCB=None,doClear=True,doDeco=True):
	vB = kv.genNorm(BMin, BMax, doLog=False, midP=None)

	if (AxCB is not None):
		AxCB.clear()
		kv.genCB(AxCB,vB,r'Radial MF B$_r$(r/r$_0)^2$ [nT]',cM=BCM,Ntk=7)
	if (doClear):
		Ax.clear()

	Br = gsph.eqNormBr(nStp)
	Ax.pcolormesh(gsph.xxi,gsph.yyi,Br,cmap=BCM,norm=vB)

	kv.SetAx(xyBds,Ax)

	if (doDeco):
		Ax.set_xlabel('X [R_S]')
		Ax.set_ylabel('Y [R_S]')
		Ax.yaxis.tick_right()
		Ax.yaxis.set_label_position('right')
	return Br


#Plot Speed at 1 AU
def PlotiSlMagV(gsph,nStp,xyBds,Ax,AxCB=None,doClear=True,doDeco=True):
	vMagV = kv.genNorm(VMin, VMax, doLog=False, midP=None)

	if (AxCB is not None):
                #Add the colorbar to AxCB
                AxCB.clear()
                kv.genCB(AxCB,vMagV,"Speed [km/s]",cM=MagVCM,Ntk=7)

	#Now do main plotting
	if (doClear):
		Ax.clear()

	V = gsph.iSliceMagV(nStp)
	lat, lon = gsph.iSliceGrid()
	Ax.pcolormesh(lon,lat,V,cmap=MagVCM,norm=vMagV)

	kv.SetAx(xyBds,Ax)

	if (doDeco):
		Ax.set_xlabel('Longitude')
		Ax.set_ylabel('Latitude')
	return V

#Plot Density at 1 AU
def PlotiSlD(gsph,nStp,xyBds,Ax,AxCB=None,doClear=True,doDeco=True):
	vD = kv.genNorm(D0Min, D0Max, doLog=False, midP=None)
	if (AxCB is not None):
		AxCB.clear()
		kv.genCB(AxCB,vD,"Number density [cm-3]",cM=D0CM,Ntk=7)

	if (doClear):
		Ax.clear()

	D = gsph.iSliceD(nStp)
	lat, lon = gsph.iSliceGrid()
	Ax.pcolormesh(lon,lat,D,cmap=D0CM,norm=vD)
	kv.SetAx(xyBds,Ax)

	if (doDeco):
		Ax.set_xlabel('Longitude')
		Ax.set_ylabel('Latitude')
		Ax.yaxis.tick_right()
		Ax.yaxis.set_label_position('right')
	return D

#Plot Br and current sheet (Br=0) at 1 AU
def PlotiSlBr(gsph,nStp,xyBds,Ax,AxCB=None,doClear=True,doDeco=True):
	vB = kv.genNorm(B0Min, B0Max, doLog=False, midP=None)
	if (AxCB is not None):
		AxCB.clear()
		kv.genCB(AxCB,vB,"Radial magnetic field [nT]",cM=BCM,Ntk=7)
	if (doClear):
		Ax.clear()

	Br = gsph.iSliceBr(nStp)
	lat, lon = gsph.iSliceGrid()
	#for contour cell-centered lon lat coordinates
	lon_c = 0.25*( lon[:-1,:-1]+lon[:-1,1:]+lon[1:,:-1]+lon[1:,1:] )
	lat_c = 0.25*( lat[:-1,:-1]+lat[:-1,1:]+lat[1:,:-1]+lat[1:,1:] )

	Ax.pcolormesh(lon,lat,Br,cmap=BCM,norm=vB)
	Ax.contour(lon_c, lat_c,Br,[0.],colors='black')
	kv.SetAx(xyBds,Ax)

	if (doDeco):
		Ax.set_xlabel('Longitude')
		Ax.set_ylabel('Latitude')
		Ax.yaxis.tick_right()
		Ax.yaxis.set_label_position('right')
		#for pic4
		Ax.set_aspect('equal')
	return Br

#Plot Br and current sheet (Br=0) at certain distance set in iSliceBr
def PlotiSlBrRotatingFrame(gsph,nStp,xyBds,Ax,AxCB=None,doClear=True,doDeco=True):
	BMin = -5.
	BMax = 5.
	vB = kv.genNorm(BMin, BMax, doLog=False, midP=None)
	if (AxCB is not None):
		AxCB.clear()
		kv.genCB(AxCB,vB,"Radial magnetic field [nT]",cM=BCM,Ntk=7)
	if (doClear):
		Ax.clear()

	#Br from the i=0
	Br = gsph.iSliceBrBound(nStp)
	lat, lon = gsph.iSliceGrid()
	
	#transform into rotating frame
	#Julian date of the initial map
	jd0 = gsph.MJDs.min()
	jd_c = gsph.MJDs[nStp]
	print (jd0, jd_c)
	#Julian date of the current solution
	time_days = (jd_c - jd0)
	print (time_days)
	omega=2*180./Tsolar

	#for contour cell-centered lon lat coordinates
	lon_c = 0.25*( lon[:-1,:-1]+lon[:-1,1:]+lon[1:,:-1]+lon[1:,1:] )
	lat_c = 0.25*( lat[:-1,:-1]+lat[:-1,1:]+lat[1:,:-1]+lat[1:,1:] )

	phi = lon_c[0,:] 
	phi_prime = (phi-omega*time_days)%(2*180.)

	if np.where(np.ediff1d(phi_prime)<0)[0].size!=0: #for the first map size =0, for other maps size=1
		ind0=np.where(np.ediff1d(phi_prime)<0)[0][0]+1
		#print 'ind = ', ind0
	else:
		ind0=0 # this is for the first map
	print('ind0 = ', ind0)

	Br = np.roll(Br, -ind0, axis = 1)

	Ax.pcolormesh(lon,lat,Br,cmap=BCM,norm=vB)
	Ax.contour(lon_c, lat_c,Br,[0.],colors='black')
	kv.SetAx(xyBds,Ax)

	if (doDeco):
		Ax.set_xlabel('Longitude')
		Ax.set_ylabel('Latitude')
		Ax.yaxis.tick_right()
		Ax.yaxis.set_label_position('right')
		#for pic4
		Ax.set_aspect('equal')
	return Br


#Plot Temperature at 1 AU
def PlotiSlTemp(gsph,nStp,xyBds,Ax,AxCB=None,doClear=True,doDeco=True):
	vT = kv.genNorm(T0Min, T0Max, doLog=False, midP=None)

	if (AxCB is not None):
		AxCB.clear()
		kv.genCB(AxCB,vT,"Temperature [MK]",cM=TCM,Ntk=7)
	if (doClear):
		Ax.clear()

	Temp = gsph.iSliceT(nStp)	
	lat, lon = gsph.iSliceGrid()
	Ax.pcolormesh(lon,lat,Temp,cmap=TCM,norm=vT)

	kv.SetAx(xyBds,Ax)

	if (doDeco):
		Ax.set_xlabel('Longitude')
		Ax.set_ylabel('Latitude')
	return Temp

#Plot Density as a function of distance
def PlotDensityProf(gsph,nStp,xyBds,Ax,AxCB=None,doClear=True,doDeco=True):
	if (doClear):
		Ax.clear()

	D = gsph.RadProfDen(nStp)
	rad  = gsph.RadialProfileGrid()

	Ax.plot(rad,D,colorProf)

	if (doDeco):
		Ax.set_xlabel('Radial distance [R_sun]')
		Ax.set_ylabel('Density [cm-3]')
		Ax.set_ylim(250.,450.)
		Ax.set_xlim(20.,220.)
                #Ax.yaxis.tick_right()
                #Ax.yaxis.set_label_position('right')
	return D

#Plot speed as a function of distance
def PlotSpeedProf(gsph,nStp,xyBds,Ax,AxCB=None,doClear=True,doDeco=True):
	if (doClear):
		Ax.clear()
	V = gsph.RadProfSpeed(nStp)
	rad  = gsph.RadialProfileGrid()
	Ax.plot(rad,V,colorProf)

	if (doDeco):
		Ax.set_xlabel('Radial distance [R_sun]')
		Ax.set_ylabel('Speed [km/s]')
		Ax.set_ylim(600.,750.)
		Ax.set_xlim(20.,220.)
	return V

def PlotFluxProf(gsph,nStp,xyBds,Ax,AxCB=None,doClear=True,doDeco=True):
	if (doClear):
		Ax.clear()
	F = gsph.RadProfFlux(nStp)
	rad  = gsph.RadialProfileGrid()
	Ax.plot(rad,F,colorProf)
	
	if (doDeco):
		Ax.set_xlabel('Radial distance [R_sun]')
		Ax.set_ylabel('RhoVr^2')
		Ax.set_ylim(180000.,280000.)
		Ax.set_xlim(20.,220.)
	return F

#Adds MPI contours
#this function is from magnetosphere Viz script. PlotMPI is not used for helio as of now 
def PlotMPI(gsph,Ax,ashd=0.5):
	gCol = mpiCol
	for i in range(gsph.Ri):
		i0 = i*gsph.dNi
		Ax.plot(gsph.xxi[i0,:],gsph.yyi[i0,:],mpiCol,linewidth=cLW,alpha=ashd)

	if (gsph.Rj>1):
		for j in range(1,gsph.Rj):
			j0 = j*gsph.dNj
			Ax.plot(gsph.xxi[:,j0], gsph.yyi[:,j0],gCol,linewidth=cLW,alpha=ashd)
			Ax.plot(gsph.xxi[:,j0],-gsph.yyi[:,j0],gCol,linewidth=cLW,alpha=ashd)
		#X-axis (+)
		Ax.plot(gsph.xxi[:,0], gsph.yyi[:,0],gCol,linewidth=cLW,alpha=ashd)
		#X-axis (-)
		j0 = (gsph.Rj)*gsph.dNj
		Ax.plot(gsph.xxi[:,j0], gsph.yyi[:,j0],gCol,linewidth=cLW,alpha=ashd)
			
