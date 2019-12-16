#!/usr/bin/env python
#Splits serial restart into MPI decomposition

import argparse
import kaipy.gamera.magsphereRescale as upscl
from argparse import RawTextHelpFormatter
import numpy as np
import h5py
import kaipy.kaiH5 as kh5

if __name__ == "__main__":
	Ri = 3
	Rj = 6
	Rk = 1	
	runid = "msphere"
	nres = "0"

	outid = "mpimsphere"

	MainS = """Splits serial Gamera restart file into MPI-decomposed restart file
	(Ri,Rj,Rk) : Output MPI decomposition
	runid/nres : Run ID string and restart number, i.e. input file = runid.Res.#nres.h5
	outid : Output Run ID
	"""

	parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
	parser.add_argument('runid',metavar='runid',help="Run ID string")
	parser.add_argument('-n',type=int,metavar="nres",default=0,help="Restart number (default: %(default)s)")
	parser.add_argument('-Ri',type=int,metavar="Ri",default=Ri,help="i-Ranks (default: %(default)s)")
	parser.add_argument('-Rj',type=int,metavar="Rj",default=Rj,help="j-Ranks (default: %(default)s)")
	parser.add_argument('-Rk',type=int,metavar="Rk",default=Rk,help="k-Ranks (default: %(default)s)")
	parser.add_argument('-o',type=str,metavar="outid",default=outid,help="Output run ID (default: %(default)s)")

	#Finalize parsing
	args = parser.parse_args()
	bStr = args.runid
	nRes = args.n
	outid = args.o
	Ri = args.Ri
	Rj = args.Rj
	Rk = args.Rk
	
	#Open input file and get data
	fIn = bStr + ".Res.%05d.h5"%(nRes)
	
	print("Reading from %s"%(fIn))
	iH5 = h5py.File(fIn,'r')
	Ns,Nv,Nk,Nj,Ni = iH5['Gas'].shape
	G = np.zeros((Ns,Nv,Nk,Nj,Ni))
	M = np.zeros((3,Nk+1,Nj+1,Ni+1))
	G[:,:,:,:,:] = iH5['Gas'][:]
	M[  :,:,:,:] = iH5['magFlux'][:]
	X = iH5['X'][:]
	Y = iH5['Y'][:]
	Z = iH5['Z'][:]

	#Now create output files
	Nkp = Nk//Rk
	Njp = Nj//Rj
	Nip = Ni//Ri

	NumG = upscl.NumG
	print("Splitting (%d,%d,%d) cells into (%d,%d,%d) x (%d,%d,%d) [Cells,MPI]"%(Ni,Nj,Nk,Nip,Njp,Nkp,Ri,Rj,Rk))
	#Loop over output slices and create restarts
	for i in range(Ri):
		for j in range(Rj):
			for k in range(Rk):
				fOut = kh5.genName(outid,i,j,k,Ri,Rj,Rk,nRes)
				
				#Open output file
				oH5 = h5py.File(fOut,'w')

				iS =  i*Nip
				iE = iS+Nip
				jS =  j*Njp
				jE = jS+Njp
				kS =  k*Nkp
				kE = kS+Nkp

				#Indices for offset ghost grid

				iSg =  iS    -NumG   + NumG #Last numG to offset to 0-based index
				iEg =  iS+Nip+NumG+1 + NumG #Last numG to offset to 0-based index
				jSg =  jS    -NumG   + NumG #Last numG to offset to 0-based index
				jEg =  jS+Njp+NumG+1 + NumG #Last numG to offset to 0-based index
				kSg =  kS    -NumG   + NumG #Last numG to offset to 0-based index
				kEg =  kS+Nkp+NumG+1 + NumG #Last numG to offset to 0-based index

				print("Writing %s"%(fOut))
				print("\tMPI (%d,%d,%d) = [%d,%d]x[%d,%d]x[%d,%d]"%(i,j,k,iS,iE,jS,jE,kS,kE))
				print("\tGrid indices = (%d,%d)x(%d,%d)x(%d,%d)"%(iSg,iEg,jSg,jEg,kSg,kEg))

				#Slice subgrids
				ijkX = X[kSg:kEg,jSg:jEg,iSg:iEg]
				ijkY = Y[kSg:kEg,jSg:jEg,iSg:iEg]
				ijkZ = Z[kSg:kEg,jSg:jEg,iSg:iEg]

				#Slice pieces out of gas and magflux
				ijkG = G[:,:,kS:kE  ,jS:jE  ,iS:iE  ]
				ijkM = M[  :,kS:kE+1,jS:jE+1,iS:iE+1]

				#Write heavy variables
				oH5.create_dataset("Gas",data=ijkG)
				oH5.create_dataset("magFlux",data=ijkM)

				#Write subgrid
				oH5.create_dataset("X",data=ijkX)
				oH5.create_dataset("Y",data=ijkY)
				oH5.create_dataset("Z",data=ijkZ)

				#Transfer attributes to output
				for ak in iH5.attrs.keys():
					aStr = str(ak)
					oH5.attrs.create(ak,iH5.attrs[aStr])

				#Close this output file
				oH5.close()
	#Close input file
	iH5.close()
