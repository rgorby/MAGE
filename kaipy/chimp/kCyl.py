#Various routines to deal with K-Cylinders from PSDs
import numpy as np
import datetime
import h5py
import kaipy.kaiViz as kv

#Get grid from K-Cyl
def getGrid(fIn,do4D=False):
	with h5py.File(fIn,'r') as hf:
		X3 = hf["X"][()].T
		Y3 = hf["Y"][()].T
		Z3 = hf["Z"][()].T
		if (do4D):
			Ai = hf["A"][()].T
		
	xx = X3[:,:,0]
	yy = Y3[:,:,0]
	Zi = Z3[0,0,:]
	Ki = 10**Zi
	Kc = 0.5*(Ki[0:-1] + Ki[1:])
	if (do4D):
		Ac = 0.5*(Ai[0:-1] + Ai[1:])
		return xx,yy,Ki,Kc,Ai,Ac
	else:
		return xx,yy,Ki,Kc

def getSlc(fIn,nStp=0,vID="jPSD",doWrap=False):
	gID = "Step#%d"%(nStp)
	with h5py.File(fIn,'r') as hf:
		V = hf[gID][vID][()].T
	if (doWrap):
		return kv.reWrap(V)
	else:
		return V

#Pressure anisotropy
#doAsym : Px/Pz-1
#!doAsym: Px/(Px+Pz)

def PIso(fIn,nStp=0,pCut=1.0e-3,doAsym=False):
	Pxy = getSlc(fIn,nStp,"Pxy")
	Pz  = getSlc(fIn,nStp,"Pz" )
	Pk = 2*Pxy+Pz
	Nx,Ny = Pz.shape
	pR = np.zeros((Nx,Ny))

	for i in range(Nx):
		for j in range(Ny):
			
			if (Pk[i,j]>pCut and Pz[i,j]>pCut):
				if (doAsym):
					pR[i,j] = Pxy[i,j]/Pz[i,j] - 1.0
				else:
					pR[i,j] = Pxy[i,j]/(Pxy[i,j]+Pz[i,j])
			else:
				if (doAsym):
					pR[i,j] = 0.0
				else:
					pR[i,j] = np.nan
	return pR

#Equatorial grids (option for wrapping for contours)
def getEQGrid(fIn,doCenter=False,doWrap=False):
	if (doWrap):
		doCenter = True
	with h5py.File(fIn,'r') as hf:
		xx = hf["X"][()].T
		yy = hf["Y"][()].T

	if (not doCenter):
		return xx,yy

	Ngi,Ngj = xx.shape
	Ni = Ngi-1; Nj = Ngj-1
	xxc = np.zeros((Ni,Nj))
	yyc = np.zeros((Ni,Nj))

	xxc = 0.25*( xx[0:Ngi-1,0:Ngj-1] + xx[1:Ngi,0:Ngj-1] + xx[0:Ngi-1,1:Ngj] + xx[1:Ngi,1:Ngj])
	yyc = 0.25*( yy[0:Ngi-1,0:Ngj-1] + yy[1:Ngi,0:Ngj-1] + yy[0:Ngi-1,1:Ngj] + yy[1:Ngi,1:Ngj])

	if (not doWrap):
		return xxc,yyc
	else:
		return kv.reWrap(xxc),kv.reWrap(yyc)	
	

