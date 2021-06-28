#!/usr/bin/env python
#Takes one MPI-decomposed restart and spits out an upscaled MPI restart

import argparse
import kaipy.embiggenUtils as upscl
from argparse import RawTextHelpFormatter
import numpy as np
import h5py
import os
import kaipy.kaiH5 as kh5

if __name__ == "__main__":
	#Input tiling
	iRi = 4
	iRj = 4
	iRk = 1
	#Output tiling
	oRi = 8
	oRj = 8
	oRk = 1

	doFast = True

	inid  = "msphere"
	outid = "msphere"

	nRes = "0"

	MainS = """Upscales and retiles a Gamera MPI resart

	(iRi,iRj,iRk) : Input MPI decomposition
	(oRi,oRj,oRk) : Output MPI decomposition
	inid/nres : Run ID string and restart number, i.e. input file = inid.MPISTUFF.Res.#nres.h5
	outid : Output Run ID
	"""	
	parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
	parser.add_argument('-i',metavar='inid',default=inid,help="Input Run ID string (default: %(default)s)")
	parser.add_argument('-n',type=int,metavar="nres",default=0,help="Restart number (default: %(default)s)")
	parser.add_argument('-iRi',type=int,metavar="iRi",default=iRi,help="Input i-Ranks (default: %(default)s)")
	parser.add_argument('-iRj',type=int,metavar="iRj",default=iRj,help="Input j-Ranks (default: %(default)s)")
	parser.add_argument('-iRk',type=int,metavar="iRk",default=iRk,help="Input k-Ranks (default: %(default)s)")
	parser.add_argument('-o',type=str,metavar="outid",default=outid,help="Output run ID (default: %(default)s)")
	parser.add_argument('-oRi',type=int,metavar="oRi",default=oRi,help="Input i-Ranks (default: %(default)s)")
	parser.add_argument('-oRj',type=int,metavar="oRj",default=oRj,help="Input j-Ranks (default: %(default)s)")
	parser.add_argument('-oRk',type=int,metavar="oRk",default=oRk,help="Input k-Ranks (default: %(default)s)")
	parser.add_argument('--down',action='store_true',default=False,help='Downscale instead of upscale (default: %(default)s)')

#Finalize parsing
	args = parser.parse_args()
	bStr = args.i
	nRes = args.n
	outid = args.o
	iRi = args.iRi
	iRj = args.iRj
	iRk = args.iRk

	oRi = args.oRi
	oRj = args.oRj
	oRk = args.oRk

	doUp = not args.down

#Start by pulling tiled restart into one brick w/ halos
	X,Y,Z,nG,nM,nB,oG,oM,oB,G0,fIn = upscl.PullRestartMPI(bStr,nRes,iRi,iRj,iRk)

#Do upscaling on all variables
	#Chop out last 2 cells on each side, then upscale

	#Start w/ grid
	rX,rY,rZ = upscl.upGrid(X,Y,Z)

	dVr = upscl.Volume(rX,rY,rZ)
	dV0 = upscl.Volume( X, Y, Z)

	#Face-centered fluxes
	nrM = upscl.upFlux(nM)

	#Now ready to do cell-centered variables
	nrG = upscl.upGas(nG,dV0,dVr,"Gas")
	nrB = upscl.upCCMag(nB,dV0,dVr,"Bxyz")
	if (G0 is not None):
		rG0 = upscl.upGas(G0,dV0,dVr,"Gas0")

	if (doFast):
		#Just replicate for oState
		orM = nrM
		orB = nrB
		orG = nrG
	else:
		orM = upscl.upFlux (oM)
		orB = upscl.upCCMag(oB,dV0,dVr,"Bxyz")
		orG = upscl.upGas  (oG,dV0,dVr,"Gas")
		
#Push back out to arbitrary decomposition
	upscl.PushRestartMPI(outid,nRes,oRi,oRj,oRk,rX,rY,rZ,nrG,nrM,nrB,orG,orM,orB,rG0,fIn)

# #Do upscaling of all variables
# 	rX = X
# 	rY = Y
# 	rZ = Z
# 	nrG = nG
# 	nrM = nM
# 	nrB = nB
# 	orG = oG
# 	orM = oM
# 	orB = oB

# 	rG0 = G0

#Toy check
	#upscl.CompRestarts(bStr,outid,nRes,iRi,iRj,iRk)
