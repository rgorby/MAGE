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
	#Close input file
	iH5.close()

	upscl.PushRestartMPI(outid,nRes,Ri,Rj,Rk,X,Y,Z,G,M,fIn)

