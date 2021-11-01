#Various routines to help with restart upscaling
#Updated to work w/ new-form restarts

import h5py
import numpy as np
import scipy
from scipy.spatial import ConvexHull
import kaipy.kaiH5 as kh5
import os
import kaipy.gamera.gamGrids as gg

TINY = 1.0e-8
NumG = 4
IDIR = 0
JDIR = 1
KDIR = 2

#Compare two restart sets to check if they're identical
def CompRestarts(iStr,oStr,nRes,Ri,Rj,Rk):
	vIDs = ["X","Y","Z","Gas","oGas","Bxyz","oBxyz","magFlux","omagFlux"]

	for i in range(Ri):
		for j in range(Rj):
			for k in range(Rk):
				fIn1 = kh5.genName(iStr,i,j,k,Ri,Rj,Rk,nRes)
				fIn2 = kh5.genName(oStr,i,j,k,Ri,Rj,Rk,nRes)

				kh5.CheckOrDie(fIn1)
				i1H5 = h5py.File(fIn1,'r')
				kh5.CheckOrDie(fIn2)
				i2H5 = h5py.File(fIn2,'r')

				#print("Rijk = %d,%d,%d"%(i,j,k))
				for v in vIDs:
					Q1 = i1H5[v][:]
					Q2 = i2H5[v][:]
					#print("\tDelta-%-12s = %12.5e, Tot = %e,%e"%(v,np.abs(Q1-Q2).sum(),np.abs(Q1).sum(),np.abs(Q2).sum()))
					

					# #Lazy interior loop
					# if (v == "magFlux"):
					# 	Nd,Nk,Nj,Ni = Q1.shape

					# 	for ii in range(Ni):
					# 		for jj in range(Nj):
					# 			for kk in range(Nk):
					# 				if (np.abs(Q1[:,kk,jj,ii]).min() > TINY):
					# 					dQ = np.abs(Q1[:,kk,jj,ii]-Q2[:,kk,jj,ii]).sum()
					# 					if (dQ>TINY):
					# 						print("ijk, dQ = %d,%d,%d,%e"%(ii,jj,kk,dQ))
					# 						print(Q1[:,kk,jj,ii])
					# 						print(Q2[:,kk,jj,ii])
									
					# 				#dQ = np.abs(Q1[:,kk,jj,ii]-Q2[:,kk,jj,ii]).sum()
					# 				#dQ = np.abs(Q1[0,kk,jj,ii]-Q2[0,kk,jj,ii]).sum()


				#Close up
				i1H5.close()
				i2H5.close()

				
#Get full data from a tiled restart file
def PullRestartMPI(bStr,nRes,Ri,Rj,Rk):
	doInit = True
	for i in range(Ri):
		for j in range(Rj):
			for k in range(Rk):
				fIn = kh5.genName(bStr,i,j,k,Ri,Rj,Rk,nRes)
				kh5.CheckOrDie(fIn)
				#print("Reading from %s"%(fIn))
				#Open input file
				iH5 = h5py.File(fIn,'r')

				if (doInit):
					fIn0 = fIn
					#Get size (w/ halos)
					Ns,Nv,Nk,Nj,Ni = iH5['Gas'].shape
					doGas0 = ('Gas0' in iH5.keys())
					Nkp = Nk-2*NumG
					Njp = Nj-2*NumG
					Nip = Ni-2*NumG

					#Get combined size (w/ ghosts)
					NkT = Rk*Nkp + 2*NumG
					NjT = Rj*Njp + 2*NumG
					NiT = Ri*Nip + 2*NumG

					#Allocate arrays (global)
					nG = np.zeros((Ns,Nv,NkT,NjT,NiT))
					oG = np.zeros((Ns,Nv,NkT,NjT,NiT))

					if (doGas0):
						G0 = np.zeros((Ns,Nv,NkT,NjT,NiT))
					else:
						G0 = None
					nM = np.zeros((3,NkT+1,NjT+1,NiT+1))
					oM = np.zeros((3,NkT+1,NjT+1,NiT+1))
					nB = np.zeros((3,NkT  ,NjT  ,NiT  ))
					oB = np.zeros((3,NkT  ,NjT  ,NiT  ))

					X = np.zeros((NkT+1,NjT+1,NiT+1))
					Y = np.zeros((NkT+1,NjT+1,NiT+1))
					Z = np.zeros((NkT+1,NjT+1,NiT+1))

					#Fill w/ nans for testing
					nM[:] = np.nan
					nG[:] = np.nan
					nB[:] = np.nan
					oM[:] = np.nan
					oG[:] = np.nan
					oB[:] = np.nan
					
					X[:] = np.nan
					Y[:] = np.nan
					Z[:] = np.nan

					doInit = False

			#Get grids
				Corner2Global(X,"X",iH5,i,j,k,Ri,Rj,Rk,Nip,Njp,Nkp)
				Corner2Global(Y,"Y",iH5,i,j,k,Ri,Rj,Rk,Nip,Njp,Nkp)
				Corner2Global(Z,"Z",iH5,i,j,k,Ri,Rj,Rk,Nip,Njp,Nkp)

			#Get cell-centered quantities
				Gas2Global(nG, "Gas",iH5,i,j,k,Ri,Rj,Rk,Nip,Njp,Nkp)
				Gas2Global(oG,"oGas",iH5,i,j,k,Ri,Rj,Rk,Nip,Njp,Nkp)
				if (doGas0):
					Gas2Global(G0,"Gas0",iH5,i,j,k,Ri,Rj,Rk,Nip,Njp,Nkp)

				Bvec2Global(nB, "Bxyz",iH5,i,j,k,Ri,Rj,Rk,Nip,Njp,Nkp)
				Bvec2Global(oB,"oBxyz",iH5,i,j,k,Ri,Rj,Rk,Nip,Njp,Nkp)
			#Get face-centered quantities
				Flux2Global(nM, "magFlux",iH5,i,j,k,Ri,Rj,Rk,Nip,Njp,Nkp)
				Flux2Global(oM,"omagFlux",iH5,i,j,k,Ri,Rj,Rk,Nip,Njp,Nkp)

				#print("\tMPI (%d,%d,%d) = [%d,%d]x[%d,%d]x[%d,%d]"%(i,j,k,iS,iE,jS,jE,kS,kE))
				#print("\tGrid indices = (%d,%d)x(%d,%d)x(%d,%d)"%(iSg,iEg,jSg,jEg,kSg,kEg))

				#Close up
				iH5.close()

# print(np.isnan(G).sum(),G.size)
# print(np.isnan(X).sum(),X.size)
# print(np.isnan(M).sum(),M.size)

	return X,Y,Z,nG,nM,nB,oG,oM,oB,G0,fIn0

def Gas2Global(G,vID,iH5,i,j,k,Ri,Rj,Rk,Nip,Njp,Nkp):
	Ns,Nv,NkT,NjT,NiT = G.shape
	Nk = Nkp + 2*NumG
	Nj = Njp + 2*NumG
	Ni = Nip + 2*NumG

	locG = np.zeros((Ns,Nv,Nk,Nj,Ni))
	locG = iH5[vID][:]

	#Set general global/local indices
	iS,iE,iSg,iEg = ccL2G(i,Ri,Nip)
	jS,jE,jSg,jEg = ccL2G(j,Rj,Njp)
	kS,kE,kSg,kEg = ccL2G(k,Rk,Nkp)

	G[:,:,kSg:kEg,jSg:jEg,iSg:iEg] = locG[:,:,kS:kE,jS:jE,iS:iE]
def Bvec2Global(B,vID,iH5,i,j,k,Ri,Rj,Rk,Nip,Njp,Nkp):
	Nd,NkT,NjT,NiT = B.shape
	Nk = Nkp + 2*NumG
	Nj = Njp + 2*NumG
	Ni = Nip + 2*NumG

	locB = np.zeros((Nd,Nk,Nj,Ni))
	locB = iH5[vID][:]

	#Set general global/local indices
	iS,iE,iSg,iEg = ccL2G(i,Ri,Nip)
	jS,jE,jSg,jEg = ccL2G(j,Rj,Njp)
	kS,kE,kSg,kEg = ccL2G(k,Rk,Nkp)

	B[:,kSg:kEg,jSg:jEg,iSg:iEg] = locB[:,kS:kE,jS:jE,iS:iE]

def Flux2Global(M,vID,iH5,i,j,k,Ri,Rj,Rk,Nip,Njp,Nkp):
	Nd,NkT1,NjT1,NiT1 = M.shape

	Nk = Nkp + 2*NumG
	Nj = Njp + 2*NumG
	Ni = Nip + 2*NumG
	locM = np.zeros((Nd,Nk+1,Nj+1,Ni+1))
	locM = iH5[vID][:]

	#Set global/local indices
	iS,iE,iSg,iEg = ccL2G(i,Ri,Nip)
	jS,jE,jSg,jEg = ccL2G(j,Rj,Njp)
	kS,kE,kSg,kEg = ccL2G(k,Rk,Nkp)
	M[:,kSg:kEg+1,jSg:jEg+1,iSg:iEg+1] = locM[:,kS:kE+1,jS:jE+1,iS:iE+1]

def ccL2G(i,Ri,Nip):
	#Set general case of only active
	iS = NumG
	iE = iS+Nip+1
	iSg = NumG + i*Nip
	iEg = iSg+Nip+1

	if (i == 0):
		#Need low-i ghosts
		iS = 0
		iSg = iSg-NumG
	if (i == (Ri-1)):
		iE  = iE +NumG
		iEg = iEg+NumG	
	return iS,iE,iSg,iEg

#Pull corner-centered variables w/ halos into global sized grid
def Corner2Global(Q,vID,iH5,i,j,k,Ri,Rj,Rk,Nip,Njp,Nkp):
	iS =  i*Nip
	iE = iS+Nip
	jS =  j*Njp
	jE = jS+Njp
	kS =  k*Nkp
	kE = kS+Nkp
	iSg =  iS    -NumG   + NumG
	iEg =  iS+Nip+NumG+1 + NumG
	jSg =  jS    -NumG   + NumG
	jEg =  jS+Njp+NumG+1 + NumG
	kSg =  kS    -NumG   + NumG
	kEg =  kS+Nkp+NumG+1 + NumG

	Q[kSg:kEg,jSg:jEg,iSg:iEg] = iH5[vID][:]

#Push restart data w/ ghosts to an output tiling
def PushRestartMPI(outid,nRes,Ri,Rj,Rk,X,Y,Z,nG,nM,nB,oG,oM,oB,G0,f0,dtScl=1.0):
	if (G0 is not None):
		doGas0 = True

	#print(X.shape,nG.shape,nM.shape,nB.shape)
	Ns,Nv,NkT,NjT,NiT = nG.shape
	Nkp = (NkT-2*NumG)//Rk
	Njp = (NjT-2*NumG)//Rj
	Nip = (NiT-2*NumG)//Ri

	print("Reading attributes from %s"%(f0))
	iH5 = h5py.File(f0,'r')

	#print("Splitting (%d,%d,%d) cells into (%d,%d,%d) x (%d,%d,%d) [Cells,MPI]"%(Ni,Nj,Nk,Nip,Njp,Nkp,Ri,Rj,Rk))
	#Loop over output slices and create restarts
	for i in range(Ri):
		for j in range(Rj):
			for k in range(Rk):
				fOut = kh5.genName(outid,i,j,k,Ri,Rj,Rk,nRes)
				#print("Writing to %s"%(fOut))

				if (os.path.exists(fOut)):
					os.remove(fOut)
				#Open output file
				oH5 = h5py.File(fOut,'w')

				#Transfer attributes to output
				for ak in iH5.attrs.keys():
					aStr = str(ak)
					#print(aStr)
					if (aStr == "dt0"):
						oH5.attrs.create(ak,dtScl*iH5.attrs[aStr])
						print(dtScl*iH5.attrs[aStr])
					else:
						
						oH5.attrs.create(ak,iH5.attrs[aStr])
			#Base indices
				iS =  i*Nip
				iE = iS+Nip
				jS =  j*Njp
				jE = jS+Njp
				kS =  k*Nkp
				kE = kS+Nkp
			#Do subgrids
				iSg =  iS    -NumG   + NumG 
				iEg =  iS+Nip+NumG+1 + NumG 
				jSg =  jS    -NumG   + NumG 
				jEg =  jS+Njp+NumG+1 + NumG 
				kSg =  kS    -NumG   + NumG 
				kEg =  kS+Nkp+NumG+1 + NumG 
				ijkX = X[kSg:kEg,jSg:jEg,iSg:iEg]
				ijkY = Y[kSg:kEg,jSg:jEg,iSg:iEg]
				ijkZ = Z[kSg:kEg,jSg:jEg,iSg:iEg]
			#Do face fluxes (same indices as subgrids)
				ijknM = nM [:,kSg:kEg,jSg:jEg,iSg:iEg]
				ijkoM = oM [:,kSg:kEg,jSg:jEg,iSg:iEg]

			#Do cell-centered values
				iSg =  iS    -NumG + NumG 
				iEg =  iS+Nip+NumG + NumG 
				jSg =  jS    -NumG + NumG 
				jEg =  jS+Njp+NumG + NumG 
				kSg =  kS    -NumG + NumG 
				kEg =  kS+Nkp+NumG + NumG

				ijknG = nG[:,:,kSg:kEg,jSg:jEg,iSg:iEg]
				ijknB = nB[:,kSg:kEg,jSg:jEg,iSg:iEg]
				ijkoG = oG[:,:,kSg:kEg,jSg:jEg,iSg:iEg]
				ijkoB = oB[:,kSg:kEg,jSg:jEg,iSg:iEg]

				if (doGas0):
					ijkG0 = G0[:,:,kSg:kEg,jSg:jEg,iSg:iEg]

			#Write vars
				oH5.create_dataset( "Gas"    ,data=ijknG)
				oH5.create_dataset( "magFlux",data=ijknM)
				oH5.create_dataset( "Bxyz"   ,data=ijknB)
				oH5.create_dataset("oGas"    ,data=ijkoG)
				oH5.create_dataset("omagFlux",data=ijkoM)
				oH5.create_dataset("oBxyz"   ,data=ijkoB)

				oH5.create_dataset("X",data=ijkX)
				oH5.create_dataset("Y",data=ijkY)
				oH5.create_dataset("Z",data=ijkZ)

				if (doGas0):
					oH5.create_dataset("Gas0",data=ijkG0)

				#Close this output file
				oH5.close()

	#Close input file
	iH5.close()

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


	return X.T,Y.T,Z.T

#Return cell centered volume from grid X,Y,Z (w/ ghosts)
def Volume(Xg,Yg,Zg):

	Ngk1,Ngj1,Ngi1 = Xg.shape
	Nk = Ngk1-1
	Nj = Ngj1-1
	Ni = Ngi1-1

	print("Calculating volume of grid of size (%d,%d,%d)"%(Ni,Nj,Nk))
	dV = np.zeros((Nk,Nj,Ni))
	ijkPts = np.zeros((8,3))

	#Assuming LFM-like symmetry
	k = 0
	for j in range(Nj):
		for i in range(Ni):
			ijkPts[:,0] = Xg[k:k+2,j:j+2,i:i+2].flatten()
			ijkPts[:,1] = Yg[k:k+2,j:j+2,i:i+2].flatten()
			ijkPts[:,2] = Zg[k:k+2,j:j+2,i:i+2].flatten()

			dV[:,j,i] = ConvexHull(ijkPts,incremental=False).volume
	return dV

#Upscale magnetic fluxes (M) on grid X,Y,Z (w/ ghosts) to doubled grid
def upFlux(M):
	
	#Chop out outer two cells, upscale inside
	cM = M[:,2:-2,2:-2,2:-2]

	Nd,Nkc,Njc,Nic = cM.shape
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
				Mu[IDIR,kp:kp+2,jp:jp+2,ip] = 0.25*cM[IDIR,k,j,i]
				#East i face (4)
				Mu[IDIR,kp:kp+2,jp:jp+2,ip+2] = 0.25*cM[IDIR,k,j,i+1]

				#South j face (4)
				Mu[JDIR,kp:kp+2,jp,ip:ip+2] = 0.25*cM[JDIR,k,j,i]
				#North j face (4)
				Mu[JDIR,kp:kp+2,jp+2,ip:ip+2] = 0.25*cM[JDIR,k,j+1,i]

				#Bottom k face (4)
				Mu[KDIR,kp,jp:jp+2,ip:ip+2] = 0.25*cM[KDIR,k,j,i]
				#Top k face (4)
				Mu[KDIR,kp+2,jp:jp+2,ip:ip+2] = 0.25*cM[KDIR,k+1,j,i]

				#Now all exterior faces are done
				#12 remaining interior faces

				#Interior i faces (4)
				Mu[IDIR,kp:kp+2,jp:jp+2,ip+1] = 0.5*( Mu[IDIR,kp:kp+2,jp:jp+2,ip] + Mu[IDIR,kp:kp+2,jp:jp+2,ip+2] )

				#Interior j faces (4)
				Mu[JDIR,kp:kp+2,jp+1,ip:ip+2] = 0.5*( Mu[JDIR,kp:kp+2,jp,ip:ip+2] + Mu[JDIR,kp:kp+2,jp+2,ip:ip+2] )

				#Interior k faces (4)
				Mu[KDIR,kp+1,jp:jp+2,ip:ip+2] = 0.5*( Mu[KDIR,kp,jp:jp+2,ip:ip+2] + Mu[KDIR,kp+2,jp:jp+2,ip:ip+2] )
	#MaxDiv(M)
	#MaxDiv(Mu)

	return Mu	
#Upscale gas variable (G) on grid X,Y,Z (w/ ghosts) to doubled grid
def upGas(G,dV0,dVu,vID="Gas"):
	#Chop out outer two cells, upscale inside
	cG = G[:,:,2:-2,2:-2,2:-2]
	cdV0 = dV0[2:-2,2:-2,2:-2]

	Ns,Nv,Nk,Nj,Ni = cG.shape
	
	Gu = np.zeros((Ns,Nv,2*Nk,2*Nj,2*Ni))
	print("Upscaling %s variables ..."%(vID))

	#Loop over coarse grid
	for s in range(Ns):
		for v in range(Nv):
			print("\tUpscaling Species %d, Variable %d"%(s,v))
			Gu[s,v,:,:,:] = upVarCC(cG[s,v,:,:,:],cdV0,dVu)

	return Gu

#Upscale Bxyz on grid X,Y,Z (w/ ghosts) to doubled grid
def upCCMag(B,dV0,dVu,vID="Bxyz"):
	#Chop out outer two cells, upscale inside
	cB = B[:,2:-2,2:-2,2:-2]
	cdV0 = dV0[2:-2,2:-2,2:-2]

	Nv,Nk,Nj,Ni = cB.shape

	Bu = np.zeros((Nv,2*Nk,2*Nj,2*Ni))
	print("Upscaling %s variables ..."%(vID))

	#Loop over coarse grid
	
	for v in range(Nv):
		print("\tUpscaling Variable %d"%(v))
		Bu[v,:,:,:] = upVarCC(cB[v,:,:,:],cdV0,dVu)

	return Bu

#Upscale single cell-centered variable Q, dV0=dV on coarse, dVu=dV on fine grid
def upVarCC(Q,dV0,dVu):

	Nk,Nj,Ni = Q.shape
	Qu = np.zeros((2*Nk,2*Nj,2*Ni))

	#Loop over coarse grid
	for k in range(Nk):
		for j in range(Nj):
			for i in range(Ni):
				QdV = Q[k,j,i]*dV0[k,j,i]
				dVijk = dVu[2*k:2*k+2,2*j:2*j+2,2*i:2*i+2] #Volumes of the finer subgrid
				vScl = dVijk.sum() #Total volume of the finer subchunks
				QdVu = QdV*(dVijk/vScl) #Give weighted contribution to each subcell
				Qu[2*k:2*k+2,2*j:2*j+2,2*i:2*i+2] = QdVu/dVijk #Scale back to density

	#Test conservation
	print("\t\tCoarse (Total) = %e"%(Q [:,:,:]*dV0).sum())
	print("\t\tFine   (Total) = %e"%(Qu[:,:,:]*dVu).sum())
	return Qu
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