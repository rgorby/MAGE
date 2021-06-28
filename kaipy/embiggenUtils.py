#Various routines to help with restart upscaling
#Updated to work w/ new-form restarts

import h5py
import numpy as np
import scipy
from scipy.spatial import ConvexHull
import kaipy.kaiH5 as kh5
import os

TINY = 1.0e-8
NumG = 4
IDIR = 0
JDIR = 1
KDIR = 2

#Get full data from a tiled restart file
def PullRestartMPI(bStr,nRes,Ri,Rj,Rk):
	doInit = True
	for i in range(Ri):
		for j in range(Rj):
			for k in range(Rk):
				fIn = kh5.genName(bStr,i,j,k,Ri,Rj,Rk,nRes)
				kh5.CheckOrDie(fIn)
				print("Reading from %s"%(fIn))
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

					#Allocate arrays
					G = np.zeros((Ns,Nv,NkT,NjT,NiT))
					if (doGas0):
						G0 = np.zeros((Ns,Nv,NkT,NjT,NiT))
					else:
						G0 = None
					M = np.zeros((3,NkT+1,NjT+1,NiT+1))
					B = np.zeros((3,NkT  ,NjT  ,NiT  ))

					X = np.zeros((NkT+1,NjT+1,NiT+1))
					Y = np.zeros((NkT+1,NjT+1,NiT+1))
					Z = np.zeros((NkT+1,NjT+1,NiT+1))
					#Fill w/ nans for testing
					M[:] = np.nan
					G[:] = np.nan
					X[:] = np.nan

					doInit = False
				#Now figure out indices of this brick
				iS =  i*Nip
				iE = iS+Nip
				jS =  j*Njp
				jE = jS+Njp
				kS =  k*Nkp
				kE = kS+Nkp

			#Do cell-centered quantities
				#Last numG to offset to 0-based index
				iSg =  iS    -NumG + NumG 
				iEg =  iS+Nip+NumG + NumG 
				jSg =  jS    -NumG + NumG 
				jEg =  jS+Njp+NumG + NumG 
				kSg =  kS    -NumG + NumG 
				kEg =  kS+Nkp+NumG + NumG

				G[:,:,kSg:kEg,jSg:jEg,iSg:iEg] = iH5['Gas'][:]
				if (doGas0):
					G0[:,:,kSg:kEg,jSg:jEg,iSg:iEg] = iH5['Gas0'][:]
				B[:,kSg:kEg,jSg:jEg,iSg:iEg] = iH5['Bxyz']

				#Do offset quantities
				iSg =  iS    -NumG   + NumG
				iEg =  iS+Nip+NumG+1 + NumG
				jSg =  jS    -NumG   + NumG
				jEg =  jS+Njp+NumG+1 + NumG
				kSg =  kS    -NumG   + NumG
				kEg =  kS+Nkp+NumG+1 + NumG

				M[:,kSg:kEg,jSg:jEg,iSg:iEg] = iH5['magFlux'][:]

				X[kSg:kEg,jSg:jEg,iSg:iEg] = iH5['X'][:]
				Y[kSg:kEg,jSg:jEg,iSg:iEg] = iH5['Y'][:]
				Z[kSg:kEg,jSg:jEg,iSg:iEg] = iH5['Z'][:]

				#print("\tMPI (%d,%d,%d) = [%d,%d]x[%d,%d]x[%d,%d]"%(i,j,k,iS,iE,jS,jE,kS,kE))
				#print("\tGrid indices = (%d,%d)x(%d,%d)x(%d,%d)"%(iSg,iEg,jSg,jEg,kSg,kEg))

				#Close up
				iH5.close()

# print(np.isnan(G).sum(),G.size)
# print(np.isnan(X).sum(),X.size)
# print(np.isnan(M).sum(),M.size)

	return X,Y,Z,G,M,B,G0,fIn0

#Push restart data w/ ghosts to an output tiling
def PushRestartMPI(outid,nRes,Ri,Rj,Rk,X,Y,Z,G,M,B,G0,f0):
	if (G0 is not None):
		doGas0 = True

	print("Reading attributes from %s"%(f0))
	iH5 = h5py.File(f0,'r')

	#print("Splitting (%d,%d,%d) cells into (%d,%d,%d) x (%d,%d,%d) [Cells,MPI]"%(Ni,Nj,Nk,Nip,Njp,Nkp,Ri,Rj,Rk))
	#Loop over output slices and create restarts
	for i in range(Ri):
		for j in range(Rj):
			for k in range(Rk):
				fOut = kh5.genName(outid,i,j,k,Ri,Rj,Rk,nRes)
				print("Writing to %s"%(fOut))

				if (os.path.exists(fOut)):
					os.remove(fOut)
				#Open output file
				oH5 = h5py.File(fOut,'w')

				#Transfer attributes to output
				for ak in iH5.attrs.keys():
					aStr = str(ak)
					oH5.attrs.create(ak,iH5.attrs[aStr])

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

			#Do cell-centered values
				iSg =  0
				iEg =  0
				jSg =  0
				jEg =  0
				kSg =  0
				kEg =  0
				ijkG  = G [:,:,kS:kE  ,jS:jE  ,iS:iE  ]
				ijkB  = B [  :,kS:kE  ,jS:jE  ,iS:iE  ]

				if (doGas0):
					ijkG0 = G0[:,:,kS:kE  ,jS:jE  ,iS:iE  ] 

			#Do face fluxes
				iSg =  0
				iEg =  0
				jSg =  0
				jEg =  0
				kSg =  0
				kEg =  0

				ijkM = M [  :,kS:kE+1,jS:jE+1,iS:iE+1]

			#Write vars
				oH5.create_dataset( "Gas"    ,data=ijkG)
				oH5.create_dataset( "magFlux",data=ijkM)
				oH5.create_dataset( "Bxyz"   ,data=ijkB)
				oH5.create_dataset("oGas"    ,data=ijkG)
				oH5.create_dataset("omagFlux",data=ijkM)
				oH5.create_dataset("oBxyz"   ,data=ijkB)

				oH5.create_dataset("X",data=ijkX)
				oH5.create_dataset("Y",data=ijkY)
				oH5.create_dataset("Z",data=ijkZ)

				if (doGas0):
					oH5.create_dataset("Gas0",data=ijkG0)

				#Close this output file
				oH5.close()

	#Close input file
	iH5.close()
