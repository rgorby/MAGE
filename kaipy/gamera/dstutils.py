import numpy as np
import kaipy.kaiH5 as kh5
import h5py

mu0 = 4*np.pi*1e-7
Mp = 1.67e-27 #[kg]
x0 = 1*6.38e6 #[m]   - RE
u0 = 100e3    #[m/s] - 100 km/s
t0 = x0/u0 #[s]	  - 
d0 = Mp*1e6 # [kg/m^3] - 1 particle/cc
p0 = d0*u0*u0 #[N/m^2]
B0 = np.sqrt(mu0*d0*u0*u0)*1e9 # [nT]


Xs = [1,0,-1,0]
Ys = [0,1,0,-1]
i0 = 1
Np = 4

def GetSteps(fIn):
        nSteps,sIds = kh5.cntSteps(fIn)
        s0 = sIds.min()
        sE = sIds.max()
        return s0,sE

def GetCC(fIn):
	X  = kh5.PullVar(fIn,"X")
	Y  = kh5.PullVar(fIn,"Y")
	Z  = kh5.PullVar(fIn,"Z")
	dV = kh5.PullVar(fIn,"dV")

	Xc,Yc,Zc = cGrid(X,Y,Z)
	return Xc,Yc,Zc,dV

def cGrid(X,Y,Z):
	Nig,Njg,Nkg = X.shape
	xyzii = np.zeros((Nig  ,Njg  ,Nkg  ,3))
	xyzcc = np.zeros((Nig-1,Njg-1,Nkg-1,3))
	xyzii[:,:,:,0] = X
	xyzii[:,:,:,1] = Y
	xyzii[:,:,:,2] = Z
	for i in range(Nig-1):
		for j in range(Njg-1):
			for k in range(Nkg-1):
				for d in range(3):
					xyzcc[i,j,k,d] = 0.125*( xyzii[i,j,k,d] + xyzii[i+1,j,k,d] \
						       + xyzii[i,j+1,k,d] + xyzii[i,j,k+1,d] \
						       + xyzii[i+1,j+1,k,d] + xyzii[i,j+1,k+1,d] \
						       + xyzii[i+1,j,k+1,d] + xyzii[i+1,j+1,k+1,d] )
	return xyzcc[:,:,:,0],xyzcc[:,:,:,1],xyzcc[:,:,:,2]


def GetBz(fIn,nSlc,Xc,Yc,Zc,dV,X=0,Y=0,Z=0):
	Rc = np.sqrt( (Xc-X)**2.0+(Yc-Y)**2.0+(Zc-Z)**2.0)
	Jx = kh5.PullVar(fIn,"Jx",nSlc)
	Jy = kh5.PullVar(fIn,"Jy",nSlc)

	Bz = -(Jx*(Yc-Y) - Jy*(Xc-X))/(Rc**3.0)
	Bz = B0*Bz*dV/(4*np.pi)
	t = kh5.tStep(fIn,nSlc)
	mjd = kh5.tStep(fIn,nSlc,"MJD")
	return mjd,t,Bz

#Get Bz from MPI data
def GetBzMPI(gsph,nSlc,Xc,Yc,Zc,dV,X=0,Y=0,Z=0):
	Rc = np.sqrt( (Xc-X)**2.0+(Yc-Y)**2.0+(Zc-Z)**2.0)
	Jx = gsph.GetVar("Jx",nSlc)
	Jy = gsph.GetVar("Jy",nSlc)
	Bz = -(Jx*Yc - Jy*Xc)/(Rc**3.0)
	Bz = B0*Bz*dV/(4*np.pi)
	t = gsph.T[nSlc-gsph.s0]
	mjd = gsph.MJDs[nSlc-gsph.s0]
	return mjd,t,Bz
	
def GetDST(fIn,nSlc,Xc,Yc,Zc,dV):
	#Get zero
	t,Bz0 = GetBz(fIn,nSlc,Xc,Yc,Zc,dV)
	d0 = Bz0[i0:,:,:].sum()
	dA = 0.0
	for n in range(Np):
		t,BzN = GetBz(fIn,nSlc,Xc,Yc,Zc,dV,Xs[n],Ys[n],0.0)
		dA = dA + BzN[i0:,:,:].sum()

	dA = dA/Np
	return t,d0,dA
def GetSymH(fBC):
	kh5.CheckOrDie(fBC)
	with h5py.File(fBC,'r') as hf:
		tData = hf['T'][()]
		dstData = hf['symh'][()]
		utData = hf['UT'][()]
	return utData,tData,dstData		
