#!/usr/bin/env python
#Generates HDF-5 grid for Gamera
#

#Grid types (gID)
#0: LFM
#1: Egg
#2: Spherical
#3: Ellipse
import argparse
import kaipy.gamera.gamGrids as gg
from argparse import RawTextHelpFormatter
import numpy as np
if __name__ == "__main__":
	#Defaults
	en = 1
	gID = 0

	doEpsY = True #Whether to set pole y values

	Rin = 2.5
	Rout = 27.5
	xtail = 250.0
	TINY = 1.0e-8 #Default tiny value
	rngtol = 1.2 #Default ring tolerance
	djwarp = -1.5

	#Rout = 29.5 #LFM to include box

	Ni0 = 48
	Nj0 = 24
	Nk0 = 32
	fOut = "Grid.h5"
	fIn = "lfm.hdf"
	MainS = """Generates HDF5 grid for Gamera
	(Ni,Nj,Nk) = en x (Ni0,Nj0,Nk0)
	Grid types (gid)
	0: LFM
	1: Egg
	2: Spherical
	3: Warp egg
	4: Fat egg"""
	#3: Ellipse"""
			
	parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
	parser.add_argument('-gid',type=int,metavar="gid",default=gID,help="Grid type ID (default: %(default)s)")
	parser.add_argument('-ni0',type=int,metavar="Ni0",default=Ni0,help="Number of i0 cells (default: %(default)s)")
	parser.add_argument('-nj0',type=int,metavar="Nj0",default=Nj0,help="Number of j0 cells (default: %(default)s)")
	parser.add_argument('-nk0',type=int,metavar="Nk0",default=Nk0,help="Number of k0 cells (default: %(default)s)")
	parser.add_argument('-en' ,type=int,metavar="en" ,default=en ,help="Cell multiplication factor (default: %(default)s)")
	parser.add_argument('-o',type=str,metavar="file",default=fOut,help="File to write output grid to (default: %(default)s)")
	parser.add_argument('-i',type=str,metavar="file",default=fIn ,help="Input LFM HDF4 file (default: %(default)s)")
	parser.add_argument('-Rin',type=float,metavar="Rin",default=Rin ,help="Inner radius (default: %(default)s)")
	parser.add_argument('-Rout',type=float,metavar="Rout",default=Rout ,help="Sunward outer radius (default: %(default)s)")
	parser.add_argument('-xtail',type=float,metavar="xtail",default=xtail ,help="Tailward outer radius (default: %(default)s)")
	parser.add_argument('-eps',type=float,metavar="eps",default=TINY ,help="Tiny number (default: %(default)s)")
	parser.add_argument('-viz', action='store_true', default=False,help="Show 2D figure of grid (default: %(default)s)")
	parser.add_argument('-chimp', action='store_true', default=False,help="Store grid in CHIMP format (default: %(default)s)")
	parser.add_argument('-rngtol',type=float,metavar="rngtol",default=rngtol ,help="Ring-avg tolerance (default: %(default)s)")
	parser.add_argument('-djwarp',type=float,metavar="djwarp",default=djwarp ,help="Theta-stretching for warped egg (default: %(default)s)")
	parser.add_argument('-v', action='store_true', default=False,help="Verbose output (default: %(default)s)")

	#Finalize parsing
	args = parser.parse_args()

	Ni0 = args.ni0
	Nj0 = args.nj0
	Nk0 = args.nk0
	en = args.en
	gID = args.gid
	fOut = args.o
	fIn = args.i
	Rin = args.Rin
	Rout = args.Rout
	TINY = args.eps
	doViz = args.viz
	doChimp = args.chimp
	rngtol = args.rngtol
	djwarp = args.djwarp
	doVerb = args.v
	xtail = np.abs(args.xtail)

	#---------------------
	#Do work
	Ni = Ni0*en
	Nj = Nj0*en
	Nk = Nk0*en

	#fO = "Grid%d.h5"%(Nk)
	
	if (gID == 0):
		print("Generating LFM grid ...")
		#print("\tReading from %s"%fIn)
		XX,YY = gg.genLFM(Ni=Ni,Nj=Nj,Rin=Rin,Rout=Rout,fIn=fIn,TINY=TINY)

	if (gID == 1):
		print("Generating Egg grid ...")
		XX,YY = gg.genEgg(Ni=Ni,Nj=Nj,Rin=Rin,Rout=Rout,xtail=xtail,TINY=TINY,A=0.0)
		
	if (gID == 2):
		print("Generating Spherical grid ...")
		#Rin = 7.5
		#Rout = 80.0
		XX,YY = gg.genSph(Ni=Ni,Nj=Nj,Rin=Rin,Rout=Rout,TINY=TINY)

	if (gID == 3):
		print("Generating Stretched Egg grid ...")
		XX,YY = gg.genEgg(Ni=Ni,Nj=Nj,Rin=Rin,Rout=Rout,xtail=xtail,TINY=TINY,A=djwarp)
		
	if (gID == 4):
		print("Generating Fat Egg grid ...")
		XX,YY = gg.genFatEgg(Ni=Ni,Nj=Nj,Rin=Rin,Rout=Rout,xtail=xtail,TINY=TINY,A=djwarp)
		

	#Calculate real outer radii, sunward/anti
	rOutS  = np.sqrt(XX[-1,0]**2.0  + YY[-1,0]**2.0)
	rOutAS = np.sqrt(XX[-1,-1]**2.0 + YY[-1,-1]**2.0)

	print("\tOutput: %s"%fOut)
	print("\tSize: (%d,%d,%d)"%(Ni,Nj,Nk))
	print("\tInner Radius: %f"%Rin)
	print("\tSunward Outer Radius: %f"%rOutS)
	print("\tTail Outer Radius: %f"%rOutAS)
	llBC = np.arcsin(np.sqrt(1.0/Rin))*180.0/np.pi

	print("\tLow-lat BC: %f"%(llBC))
		
	print("\nWriting to %s"%(fOut))

	xxG,yyG = gg.Aug2D(XX,YY,doEps=doEpsY)
	X3,Y3,Z3 = gg.Aug3D(xxG,yyG,Nk=Nk,TINY=TINY)

	#Debugging stuff
	# rr0 = np.sqrt(xxG**2.0 + yyG**2.0)
	# pp0 = np.arctan2(yyG,xxG)
	# print(pp0[1,:])
	# print(pp0[8,:])
	
	if (doChimp):
		gg.WriteChimp(X3,Y3,Z3,fOut=fOut)
	else:
		#Do ring checking
		gg.genRing(XX,YY,Nk=Nk,Tol=rngtol,doVerb=doVerb)
		gg.WriteGrid(X3,Y3,Z3,fOut=fOut)
	if (doViz):
		gg.VizGrid(XX,YY,xxG,yyG,fOut=fOut)
