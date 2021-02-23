#Various helper routines for heliosphere quick look plots

import argparse
from argparse import RawTextHelpFormatter
import matplotlib as mpl
import numpy as np

import kaipy.kaiViz as kv
import kaipy.gamhelio.heliosphere as hsph
import os

VMax = 800.
VMin = 300.
MagVCM = "inferno"

DMax = 150.
DMin = 2000.
DCM = "copper_r"

D0Max = 15.
D0Min = 1.
D0CM = "copper_r"

TMax = 0.2
TMin = 0.01
TCM = "copper"

BMax = 150.
BMin = -150.
BCM = "coolwarm"

dbMax = 25.0
dbCM = "RdGy_r"
bzCM = "bwr"
cLW = 0.25
bz0Col = "magenta"
mpiCol = "deepskyblue"

jMax = 10.0 #Max current for contours

#Default pressure colorbar
vP = kv.genNorm(1.0e-2,10.0,doLog=True)

#Add different size options to argument

#not used for helio right now
def AddSizeArgs(parser):
	parser.add_argument('-small' , action='store_true', default=False,help="Use smaller domain bounds (default: %(default)s)")
	parser.add_argument('-big'   , action='store_true', default=False,help="Use larger domain bounds (default: %(default)s)")
	parser.add_argument('-bigger', action='store_true', default=False,help="Use larger-er domain bounds (default: %(default)s)")
	parser.add_argument('-huge'  , action='store_true', default=False,help="Use huge domain bounds (default: %(default)s)")

#Return domain size from parsed arguments; see msphViz for options
def GetSizeBds():
	xMax = 216.
	yMax = 216.
	xyBds = [-xMax,xMax,-yMax,yMax]

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

	#Ax.contour(kv.reWrap(gsph.xxc),kv.reWrap(gsph.yyc),kv.reWrap(Bz),[0.0],colors=bz0Col,linewidths=cLW)
	kv.SetAx(xyBds,Ax)

	if (doDeco):
		#kv.addEarth2D(ax=Ax)
		Ax.set_xlabel('X [R_S]')
		Ax.set_ylabel('Y [R_S]')
	return MagV

def PlotMerMagV(gsph,nStp,xyBds,Ax,AxCB=None,doClear=True,doDeco=True):
        vMagV = kv.genNorm(VMin, VMax, doLog=False, midP=None)

        if (AxCB is not None):
                #Add the colorbar to AxCB
                AxCB.clear()
                kv.genCB(AxCB,vMagV,"Speed [km/s]",cM=MagVCM,Ntk=7)

        if (doClear):
                Ax.clear()

        Vr, Vl = gsph.MerMagV(nStp)
        xr, zr, xl, zl, r = gsph.MeridGridHalfs()
        Ax.pcolormesh(xr,zr,Vr,cmap=MagVCM,norm=vMagV)
        Ax.pcolormesh(xl,zl,Vl,cmap=MagVCM,norm=vMagV)

        kv.SetAx(xyBds,Ax)

        if (doDeco):
                Ax.set_xlabel('X [R_S]')
                Ax.set_ylabel('Z [R_S]')
        return Vr, Vl


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

def PlotMerBrNorm(gsph,nStp,xyBds,Ax,AxCB=None,doClear=True,doDeco=True):
	vB = kv.genNorm(BMin, BMax, doLog=False, midP=None)

	if (AxCB is not None):
		#Add the colorbar to AxCB
		AxCB.clear()
		kv.genCB(AxCB,vB,r'Radial MF $B_r(r/r_0)^2$ [nT]',cM=BCM,Ntk=7)
	if (doClear):
		Ax.clear()

	Br_r, Br_l = gsph.MerBrNrm(nStp)
	xr, zr, xl, zl, r = gsph.MeridGridHalfs()
	Ax.pcolormesh(xr,zr,Br_r,cmap=BCM,norm=vB,shading='auto')
	Ax.pcolormesh(xl,zl,Br_l,cmap=BCM,norm=vB,shading='auto')

	kv.SetAx(xyBds,Ax)

	if (doDeco):
		Ax.set_xlabel('X [R_S]')
		Ax.set_ylabel('Z [R_S]')
		Ax.yaxis.tick_right()
		Ax.yaxis.set_label_position('right')
	return Br_r, Br_l

def PlotMerTemp(gsph,nStp,xyBds,Ax,AxCB=None,doClear=True,doDeco=True):
	vT = kv.genNorm(TMin, TMax, doLog=False, midP=None)

	if (AxCB is not None):
		AxCB.clear()
		kv.genCB(AxCB,vT,"Temperature [MK]",cM=TCM,Ntk=7)
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

#Plot normalized density in equatorial plane
def PlotEqD(gsph,nStp,xyBds,Ax,AxCB=None,doClear=True,doDeco=True):
	vD = kv.genNorm(DMin, DMax, doLog=False, midP=None)
	
	if (AxCB is not None):
		#Add the colorbar to AxCB
		AxCB.clear()
		kv.genCB(AxCB,vD,r"Density n$(r/r_0)^2$ [cm$^{-3}$]",cM=DCM,Ntk=7)

	#Now do main plotting
	if (doClear):
		Ax.clear()

	NormD = gsph.eqNormD(nStp)
	Ax.pcolormesh(gsph.xxi,gsph.yyi,NormD,cmap=DCM,norm=vD)

	#Ax.contour(kv.reWrap(gsph.xxc),kv.reWrap(gsph.yyc),kv.reWrap(Bz),[0.0],colors=bz0Col,linewidths=cLW)
	kv.SetAx(xyBds,Ax)

	if (doDeco):
		#kv.addEarth2D(ax=Ax)
		Ax.set_xlabel('X [R_S]')
		Ax.set_ylabel('Y [R_S]')
		Ax.yaxis.tick_right()
		Ax.yaxis.set_label_position('right')
	return NormD

#Plot Temperature in equatorial plane
def PlotEqTemp(gsph,nStp,xyBds,Ax,AxCB=None,doClear=True,doDeco=True):
	vT = kv.genNorm(TMin, TMax, doLog=False, midP=None)

	if (AxCB is not None):
		AxCB.clear()
		kv.genCB(AxCB,vT,"Temperature [MK]",cM=TCM,Ntk=7)
	if (doClear):
		Ax.clear()

	Temp = gsph.eqTemp(nStp)
	Ax.pcolormesh(gsph.xxi,gsph.yyi,Temp,cmap=TCM,norm=vT)
	
	kv.SetAx(xyBds,Ax)

	if (doDeco):
		Ax.set_xlabel('X [R_S]')
		Ax.set_ylabel('Y [R_S]')
	return Temp

#Plor Br in eq plane
def PlotEqBr(gsph,nStp,xyBds,Ax,AxCB=None,doClear=True,doDeco=True):
	vB = kv.genNorm(BMin, BMax, doLog=False, midP=None)

	if (AxCB is not None):
		AxCB.clear()
		kv.genCB(AxCB,vB,r'Radial MF $B_r(r/r_0)^2$ [nT]',cM=BCM,Ntk=7)
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



def PlotiSlMagV(gsph,nStp,xyBds,Ax,AxCB=None,doClear=True,doDeco=True):
	xyBds = [0.,360.,-75.,75.]
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

        #Ax.contour(kv.reWrap(gsph.xxc),kv.reWrap(gsph.yyc),kv.reWrap(Bz),[0.0],colors=bz0Col,linewidths=cLW)
	kv.SetAx(xyBds,Ax)

	if (doDeco):
		Ax.set_xlabel('Longitude')
		Ax.set_ylabel('Latitude')
	return V

#plot density at 1 AU (last cell)
def PlotiSlD(gsph,nStp,xyBds,Ax,AxCB=None,doClear=True,doDeco=True):
	xyBds = [0.,360.,-75.,75.]
	vD = kv.genNorm(D0Min, D0Max, doLog=False, midP=None)
	if (AxCB is not None):
		AxCB.clear()
		kv.genCB(AxCB,vD,"Number density [cm-3]",cM=D0CM,Ntk=7)

	#Now do main plotting
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

def PlotiSlBr(gsph,nStp,xyBds,Ax,AxCB=None,doClear=True,doDeco=True):
	xyBds = [0.,360.,-75.,75.]
	BMin = -5.
	BMax = 5.
	vB = kv.genNorm(BMin, BMax, doLog=False, midP=None)
	if (AxCB is not None):
		AxCB.clear()
		kv.genCB(AxCB,vB,"Radial magnetic field [nT]",cM=BCM,Ntk=7)
	if (doClear):
		Ax.clear()

	Br = gsph.iSliceBr(nStp)
	lat, lon = gsph.iSliceGrid()
	Ax.pcolormesh(lon,lat,Br,cmap=BCM,norm=vB)
	kv.SetAx(xyBds,Ax)

	if (doDeco):
		Ax.set_xlabel('Longitude')
		Ax.set_ylabel('Latitude')
		Ax.yaxis.tick_right()
		Ax.yaxis.set_label_position('right')
	return Br

def PlotiSlTemp(gsph,nStp,xyBds,Ax,AxCB=None,doClear=True,doDeco=True):
	xyBds = [0.,360.,-75.,75.]
	TMin = 0.02
	TMax = 0.12
	vT = kv.genNorm(TMin, TMax, doLog=False, midP=None)

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


#Plot equatorial field
def PlotEqB(gsph,nStp,xyBds,Ax,AxCB=None,doClear=True,doDeco=True,doBz=False):
	vBZ = kv.genNorm(dbMax)
	vDB = kv.genNorm(dbMax)
	
	if (AxCB is not None):
		#Add the colorbar to AxCB
		AxCB.clear()
		if (doBz):
			kv.genCB(AxCB,vBZ,"Vertical Field [nT]",cM=bzCM,Ntk=7)
		else:	
			kv.genCB(AxCB,vDB,"Residual Field [nT]",cM=dbCM,Ntk=7)
	#Now do main plotting
	if (doClear):
		Ax.clear()
	Bz = gsph.EggSlice("Bz",nStp,doEq=True)
	if (doBz):
		Ax.pcolormesh(gsph.xxi,gsph.yyi,Bz,cmap=bzCM,norm=vBZ)
	else:	
		dbz = gsph.DelBz(nStp)
		Ax.pcolormesh(gsph.xxi,gsph.yyi,dbz,cmap=dbCM,norm=vDB)
	Ax.contour(kv.reWrap(gsph.xxc),kv.reWrap(gsph.yyc),kv.reWrap(Bz),[0.0],colors=bz0Col,linewidths=cLW)

	kv.SetAx(xyBds,Ax)

	if (doDeco):
		kv.addEarth2D(ax=Ax)
		Ax.set_xlabel('SM-X [Re]')
		Ax.set_ylabel('SM-Y [Re]')
	return Bz
	
def PlotMerid(gsph,nStp,xyBds,Ax,doDen=False,doRCM=False,AxCB=None,doClear=True,doDeco=True):
	CMx = "viridis"
	if (doDen):
		cbStr = "Density [#/cc]"
		if (doRCM):
			vN = kv.genNorm(1.0,1.0e+3,doLog=True)
		else:
			vN = kv.genNorm(0,25)
		Q = gsph.EggSlice("D",nStp,doEq=False)
	else:
		cbStr = "Pressure [nPa]"
		vN = vP
		Q = gsph.EggSlice("P",nStp,doEq=False)
	if (AxCB is not None):
		#Add the colorbar to AxCB
		AxCB.clear()
		kv.genCB(AxCB,vN,cbStr,cM=CMx)
	Ax.pcolormesh(gsph.xxi,gsph.yyi,Q,cmap=CMx,norm=vN)

	kv.SetAx(xyBds,Ax)
	if (doDeco):
		kv.addEarth2D(ax=Ax)
		Ax.set_xlabel('SM-X [Re]')
		Ax.set_ylabel('SM-Z [Re]')
		Ax.yaxis.tick_right()
		Ax.yaxis.set_label_position('right')

def PlotJyXZ(gsph,nStp,xyBds,Ax,AxCB=None,jScl=None,doDeco=True):
	if (jScl is None):
		#VERY LAZY scaling for current
		#Scale current to nA/m^2
		#Current is curl(B) = mag field / Re
		# mu0 J = curl(B)
		jScl = 4.58/( (6.371*1.0e+6)*(4.0*np.pi*1.0e+2) )  # Convert field to nT/m
		#jScl => A/m2
		jScl = jScl*1.0e+9
	vJ = kv.genNorm(jMax)
	jCMap = "PRGn"
	Nc = 15
	cVals = np.linspace(-jMax,jMax,Nc)

	if (AxCB is not None):
		AxCB.clear()
		kv.genCB(AxCB,vJ,"Jy [nA/m2]",cM=jCMap)
	Q = jScl*gsph.EggSlice("Jy",nStp,doEq=False)
	#Zero out first shell b/c bad derivative
	print(Q.shape)
	Q[0:2,:] = 0.0
	#Ax.contour(kv.reWrap(gsph.xxc),kv.reWrap(gsph.yyc),kv.reWrap(Q),cVals,norm=vJ,cmap=jCMap,linewidths=cLW)
	Ax.pcolormesh(gsph.xxi,gsph.yyi,Q,norm=vJ,cmap=jCMap)
	kv.SetAx(xyBds,Ax)
	if (doDeco):
		kv.addEarth2D(ax=Ax)
		Ax.set_xlabel('SM-X [Re]')
		Ax.set_ylabel('SM-Z [Re]')
		Ax.yaxis.tick_right()
		Ax.yaxis.set_label_position('right')

#Add MPI contours
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
			
