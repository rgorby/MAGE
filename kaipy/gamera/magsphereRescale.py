#Various routines to help with restart upscaling

import h5py
import numpy as np
import kaipy.gamera.gamGrids as gg
import scipy
from scipy.spatial import ConvexHull

TINY = 1.0e-8
NumG = 4
IDIR = 0
JDIR = 1
KDIR = 2

#Generate name of restart file
def genName(bStr,i,j,k,Ri,Rj,Rk,nRes):
	n = j + i*Rj + k*Ri*Rj
	fID = bStr + "_%04d_%04d_%04d_%04d_%04d_%04d_%012d"%(Ri,Rj,Rk,i,j,k,n)+".Res.%05d.h5"%(nRes)
	return fID

#Downscale a grid (with ghosts, k-j-i order)
def downGrid(X,Y,Z):
	Ngk,Ngj,Ngi = X.shape
	
	Nk = Ngk-2*NumG-1
	Nj = Ngj-2*NumG-1
	Ni = Ngi-2*NumG-1

	#Assuming LFM-style grid, get upper half plane
	
	xx = X[NumG,NumG:-NumG,NumG:-NumG].T
	yy = Y[NumG,NumG:-NumG,NumG:-NumG].T

	#Now half cells in r-phi polar
	rr = np.sqrt(xx**2.0+yy**2.0)
	pp = np.arctan2(yy,xx)

	NiD = Ni/2
	NjD = Nj/2
	NkD = Nk/2

	rrD = np.zeros((NiD+1,NjD+1))
	ppD = np.zeros((NiD+1,NjD+1))

	for i in range(NiD+1):
		for j in range(NjD+1):
			rrD[i,j] = rr[2*i,2*j]
			ppD[i,j] = pp[2*i,2*j]

	#Convert back to 2D Cartesian
	xxD = rrD*np.cos(ppD)
	yyD = rrD*np.sin(ppD)

	#Augment w/ 2D ghosts
	xxG,yyG = gg.Aug2D(xxD,yyD,doEps=True,TINY=TINY)

	#Augment w/ 3D ghosts
	X,Y,Z = gg.Aug3D(xxG,yyG,Nk=NkD,TINY=TINY)

	return X,Y,Z
#Upscale a grid (with ghosts, k-j-i order)
def upGrid(X,Y,Z):
	Ngk,Ngj,Ngi = X.shape
	
	Nk = Ngk-2*NumG-1
	Nj = Ngj-2*NumG-1
	Ni = Ngi-2*NumG-1
	

	#Assuming LFM-style grid, get upper half plane
	
	xx = X[NumG,NumG:-NumG,NumG:-NumG].T
	yy = Y[NumG,NumG:-NumG,NumG:-NumG].T

	#Now half cells in r-phi polar
	rr = np.sqrt(xx**2.0+yy**2.0)
	pp = np.arctan2(yy,xx)

	NiH = 2*Ni
	NjH = 2*Nj
	NkH = 2*Nk

	rrH = np.zeros((NiH+1,NjH+1))
	ppH = np.zeros((NiH+1,NjH+1))

	#Embed old points into new grid
	for i in range(Ni+1):
		for j in range(Nj+1):
			rrH[2*i,2*j] = rr[i,j]
			ppH[2*i,2*j] = pp[i,j]

	#Create I midpoints
	for i in range(Ni):
		rrH[2*i+1,:] = 0.5*( rrH[2*i,:] + rrH[2*i+2,:] )
		ppH[2*i+1,:] = 0.5*( ppH[2*i,:] + ppH[2*i+2,:] )

	#Create J midpoints
	for j in range(Nj):
		rrH[:,2*j+1] = 0.5*( rrH[:,2*j] + rrH[:,2*j+2])
		ppH[:,2*j+1] = 0.5*( ppH[:,2*j] + ppH[:,2*j+2])

	#Create I-J midpoints
	for i in range(Ni):
		for j in range(Nj):
			rrH[2*i+1,2*j+1] = 0.25*( rrH[2*i,2*j] + rrH[2*i,2*j+2] + rrH[2*i+2,2*j] + rrH[2*i+2,2*j+2] )
			ppH[2*i+1,2*j+1] = 0.25*( ppH[2*i,2*j] + ppH[2*i,2*j+2] + ppH[2*i+2,2*j] + ppH[2*i+2,2*j+2] )


	#Convert back to 2D Cartesian
	xxH = rrH*np.cos(ppH)
	yyH = rrH*np.sin(ppH)

	#Augment w/ 2D ghosts
	xxG,yyG = gg.Aug2D(xxH,yyH,doEps=True,TINY=TINY)

	#Augment w/ 3D ghosts
	X,Y,Z = gg.Aug3D(xxG,yyG,Nk=NkH,TINY=TINY)


	return X,Y,Z

#Downscale gas variable (G) on grid X,Y,Z (w/ ghosts) to halved grid
def downGas(X,Y,Z,G,Xd,Yd,Zd):
	Ns,Nv,Nk,Nj,Ni = G.shape
	dV  = Volume(X ,Y ,Z )
	dVd = Volume(Xd,Yd,Zd)
	print("Volume ratio (Coarse/Fine) = %f"%(dVd.sum()/dV.sum()))
	Gd = np.zeros((Ns,Nv,Nk/2,Nj/2,Ni/2))
	print("Downscaling gas variables ...")

	#Loop over coarse grid
	for s in range(Ns):
		for v in range(Nv):
			print("\tDownscaling Species %d, Variable %d"%(s,v))
			for k in range(Nk/2):
				for j in range(Nj/2):
					for i in range(Ni/2):
						dVijk = dV[2*k:2*k+2,2*j:2*j+2,2*i:2*i+2] #Volumes of the finer subgrid
						dQijk = G[s,v,2*k:2*k+2,2*j:2*j+2,2*i:2*i+2] #stuff density in finer subgrid
						#Stuff in the subchunk
						QdV = (dVijk*dQijk).sum()
						Gd[s,v,k,j,i] = QdV/dVd[k,j,i] #Scale back to density
			#Test conservation
			print("\t\tCoarse (Total) = %e"%(Gd[s,v,:,:,:]*dVd).sum())
			print("\t\tFine   (Total) = %e"%(G [s,v,:,:,:]*dV ).sum())

	return Gd
#Upscale gas variable (G) on grid X,Y,Z (w/ ghosts) to doubled grid
def upGas(X,Y,Z,G,Xu,Yu,Zu):
	Ns,Nv,Nk,Nj,Ni = G.shape

	dV  = Volume(X ,Y ,Z )
	dVu = Volume(Xu,Yu,Zu)
	
	print("Volume ratio (Coarse/Fine) = %f"%(dV.sum()/dVu.sum()))
	Gu = np.zeros((Ns,Nv,2*Nk,2*Nj,2*Ni))
	print("Upscaling gas variables ...")
	#Loop over coarse grid
	for s in range(Ns):
		for v in range(Nv):
			print("\tUpscaling Species %d, Variable %d"%(s,v))
			for k in range(Nk):
				for j in range(Nj):
					for i in range(Ni):
						QdV = G[s,v,k,j,i]*dV[k,j,i]
						dVijk = dVu[2*k:2*k+2,2*j:2*j+2,2*i:2*i+2] #Volumes of the finer subgrid
						vScl = dVijk.sum() #Total volume of the finer subchunks
						QdVu = QdV*(dVijk/vScl) #Give weighted contribution to each subcell
						Gu[s,v,2*k:2*k+2,2*j:2*j+2,2*i:2*i+2] = QdVu/dVijk #Scale back to density
			#Test conservation
			print("\t\tCoarse (Total) = %e"%(G [s,v,:,:,:]* dV).sum())
			print("\t\tFine   (Total) = %e"%(Gu[s,v,:,:,:]*dVu).sum())
	return Gu

#Return cell centered volume (active only) from grid X,Y,Z (w/ ghosts)
def Volume(Xg,Yg,Zg):
	Ngk,Ngj,Ngi = Xg.shape
	
	Nk = Ngk-2*NumG-1
	Nj = Ngj-2*NumG-1
	Ni = Ngi-2*NumG-1

	X = Xg[NumG:-NumG,NumG:-NumG,NumG:-NumG]
	Y = Yg[NumG:-NumG,NumG:-NumG,NumG:-NumG]
	Z = Zg[NumG:-NumG,NumG:-NumG,NumG:-NumG]

	print("Calculating volume of grid of size (%d,%d,%d)"%(Ni,Nj,Nk))


	dV = np.zeros((Nk,Nj,Ni))
	ijkPts = np.zeros((8,3))

	#Assuming LFM-like symmetry
	k = 0
	for j in range(Nj):
		for i in range(Ni):
			ijkPts[:,0] = X[k:k+2,j:j+2,i:i+2].flatten()
			ijkPts[:,1] = Y[k:k+2,j:j+2,i:i+2].flatten()
			ijkPts[:,2] = Z[k:k+2,j:j+2,i:i+2].flatten()

			dV[:,j,i] = ConvexHull(ijkPts,incremental=False).volume

	return dV
#Downscale magnetic fluxes (M) on grid X,Y,Z (w/ ghosts) to halved grid
def downFlux(X,Y,Z,M,Xu,Yu,Zu):
	Nd,Nkc,Njc,Nic = M.shape
	Nk = Nkc-1
	Nj = Njc-1
	Ni = Nic-1

	Md = np.zeros((Nd,Nk/2+1,Nj/2+1,Ni/2+1))

	#Loop over coarse grid cells
	print("Downscaling face fluxes ...")
	for k in range(Nk/2):
		for j in range(Nj/2):
			for i in range(Ni/2):
				ip = 2*i
				jp = 2*j
				kp = 2*k

				#West/east i faces
				Md[IDIR,k  ,j  ,i  ] = M[IDIR,kp:kp+2,jp:jp+2,ip  ].sum()
				Md[IDIR,k  ,j  ,i+1] = M[IDIR,kp:kp+2,jp:jp+2,ip+2].sum()

				#South/north j faces
				Md[JDIR,k  ,j  ,i  ] = M[JDIR,kp:kp+2,jp  ,ip:ip+2].sum()
				Md[JDIR,k  ,j+1,i  ] = M[JDIR,kp:kp+2,jp+2,ip:ip+2].sum()
				
				#Bottom/top k faces
				Md[KDIR,k  ,j  ,i  ] = M[KDIR,kp  ,jp:jp+2,ip:ip+2].sum()
				Md[KDIR,k+1,j  ,i  ] = M[KDIR,kp+2,jp:jp+2,ip:ip+2].sum()

	return Md

#Upscale magnetic fluxes (M) on grid X,Y,Z (w/ ghosts) to doubled grid
def upFlux(X,Y,Z,M,Xu,Yu,Zu):
	Nd,Nkc,Njc,Nic = M.shape
	Nk = Nkc-1
	Nj = Njc-1
	Ni = Nic-1

	Mu = np.zeros((Nd,2*Nk+1,2*Nj+1,2*Ni+1))

	#Loop over coarse grid
	print("Upscaling face fluxes ...")
	for k in range(Nk):
		for j in range(Nj):
			for i in range(Ni):
				ip = 2*i
				jp = 2*j
				kp = 2*k

				#West i face (4)
				Mu[IDIR,kp:kp+2,jp:jp+2,ip] = 0.25*M[IDIR,k,j,i]
				#East i face (4)
				Mu[IDIR,kp:kp+2,jp:jp+2,ip+2] = 0.25*M[IDIR,k,j,i+1]

				#South j face (4)
				Mu[JDIR,kp:kp+2,jp,ip:ip+2] = 0.25*M[JDIR,k,j,i]
				#North j face (4)
				Mu[JDIR,kp:kp+2,jp+2,ip:ip+2] = 0.25*M[JDIR,k,j+1,i]

				#Bottom k face (4)
				Mu[KDIR,kp,jp:jp+2,ip:ip+2] = 0.25*M[KDIR,k,j,i]
				#Top k face (4)
				Mu[KDIR,kp+2,jp:jp+2,ip:ip+2] = 0.25*M[KDIR,k+1,j,i]

				#Now all exterior faces are done
				#12 remaining interior faces

				#Interior i faces (4)
				Mu[IDIR,kp:kp+2,jp:jp+2,ip+1] = 0.5*( Mu[IDIR,kp:kp+2,jp:jp+2,ip] + Mu[IDIR,kp:kp+2,jp:jp+2,ip+2] )

				#Interior j faces (4)
				Mu[JDIR,kp:kp+2,jp+1,ip:ip+2] = 0.5*( Mu[JDIR,kp:kp+2,jp,ip:ip+2] + Mu[JDIR,kp:kp+2,jp+2,ip:ip+2] )

				#Interior k faces (4)
				Mu[KDIR,kp+1,jp:jp+2,ip:ip+2] = 0.5*( Mu[KDIR,kp,jp:jp+2,ip:ip+2] + Mu[KDIR,kp+2,jp:jp+2,ip:ip+2] )
	#MaxDiv(M)

	return Mu

#Calculates maximum divergence of mag flux data
def MaxDiv(M):
	Nd,Nkc,Njc,Nic = M.shape
	Nk = Nkc-1
	Nj = Njc-1
	Ni = Nic-1

	Div = np.zeros((Nk,Nj,Ni))
	for k in range(Nk):
		for j in range(Nj):
			for i in range(Ni):
				Div[k,j,i] = M[IDIR,k,j,i+1]-M[IDIR,k,j,i]+M[JDIR,k,j+1,i]-M[JDIR,k,j,i]+M[KDIR,k+1,j,i]-M[KDIR,k,j,i]

	mDiv = np.abs(Div).max()
	bDiv = np.abs(Div).mean()
	print("Max/Mean divergence = %e,%e"%(mDiv,bDiv))
	print("Sum divergence = %e"%(Div.sum()))

	return Div