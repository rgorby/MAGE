#!/usr/bin/env python
#Generates LFM-style HDF-5 grid for Gamera

import argparse
import kaipy.gamera.gamGrids as gg
from argparse import RawTextHelpFormatter
import numpy as np

#Ring params
rParams = {
	"D": '<ring gid="lfm" doRing="T" Nr="4" Nc1="8" Nc2="16" Nc3="32" Nc4="32"/>',
	"Q": '<ring gid="lfm" doRing="T" Nr="8" Nc1="8" Nc2="16" Nc3="32" Nc4="32" Nc5="64" Nc6="64" Nc7="64" Nc8="64"/>',
	"O": '<ring gid="lfm" doRing="T" Nr="12" Nc1="8" Nc2="16" Nc3="32" Nc4="32" Nc5="64" Nc6="64" Nc7="64" Nc8="64" Nc9="128" Nc10="128" Nc11="128" Nc12="128"/>',
	"H": '<ring gid="lfm" doRing="T" Nr="10" Nc1="16" Nc2="32" Nc3="64" Nc4="64" Nc5="128" Nc6="128" Nc7="128" Nc8="256" Nc9="256" Nc10="256"/>'

}

#"H": '<ring gid="lfm" doRing="T" Nr="16" Nc1="8" Nc2="16" Nc3="32" Nc4="32" Nc5="64" Nc6="64" Nc7="64" Nc8="64" Nc9="128" Nc10="128" Nc11="128" Nc12="128" Nc13="256" Nc14="256" Nc15="256" Nc16="256"/>'


if __name__ == "__main__":

#Arg parsing
	Nc0 = 8 #Number of outer i cells to cut out from LFM grid (OCT)
	fIn = "./lfmG"
	doEpsY = True
	TINY = 1.0e-8
	Rin = 2.0
	Rout = 0.0
	#List of grids
	gStrs = ['D','Q','O','H']
	gLabs = ["Double","Quad","Oct","Hex"]
	Nij0 = 48
	Nk0  = 64

	MainS = """Generates LFM-style HDF5 grid for Gamera
	Grid types (gid)
	D: Double ( 48, 48, 64)
	Q: Quad   ( 96, 96,128)
	O: Oct    (192,192,256)
	H: Hex    (384,384,512)
	"""
	
	parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
	parser.add_argument('-gid',type=str,default="D",choices=gStrs,help="Grid Resolution Specifier (default: %(default)s)")
	parser.add_argument('-viz', action='store_true', default=False,help="Show 2D figure of grid (default: %(default)s)")
	parser.add_argument('-chimp', action='store_true', default=False,help="Store grid in CHIMP format (default: %(default)s)")
	parser.add_argument('-Rin',type=float,metavar="Rin",default=Rin ,help="Inner radius (default: %(default)s)")
	parser.add_argument('-Rout',type=float,metavar="Rout",default=Rout ,help="Sunward outer radius (default: %(default)s)")
	parser.add_argument('-vizG', action='store_true', default=False,help="Show 2D figure w/ ghosts (default: %(default)s)")

	#Finalize parsing
	args = parser.parse_args()
	gid = args.gid
	doChimp = args.chimp
	doViz = args.viz
	doVizG = args.vizG
	Rin = args.Rin
	Rout = args.Rout
	if (doVizG):
		doViz = True

	n0 = gStrs.index(gid)
	en = 2**(n0)
	Nij = Nij0*en
	Nk  = Nk0 *en

	Ni = Nij
	Nj = Nij
	print("Generating %s LFM-style grid ...\n"%(gLabs[n0]))
	fOut = "lfm%s.h5"%(gStrs[n0])

	#Read tab data of LFM grid
	xx0,yy0 = gg.LoadTabG(fIn,Nc0)

	#Regrid to new dimensions
	XX,YY = gg.regrid(xx0,yy0,Nij,Nij,Rin=Rin,Rout=Rout)

	Rin = XX[0,0]
	#Calculate real outer radii, sunward/anti
	rOutS  = np.sqrt(XX[-1,0]**2.0  + YY[-1,0]**2.0)
	rOutAS = np.sqrt(XX[-1,-1]**2.0 + YY[-1,-1]**2.0)
	llBC = np.arcsin(np.sqrt(1.0/Rin))*180.0/np.pi

	#Do full grid
	xxG,yyG = gg.Aug2D(XX,YY,doEps=doEpsY,TINY=TINY)
	X3,Y3,Z3 = gg.Aug3D(xxG,yyG,Nk=Nk,TINY=TINY)

	#Write grid
	if (doChimp):
		gg.WriteChimp(X3,Y3,Z3,fOut=fOut)
	else:
		gg.WriteGrid(X3,Y3,Z3,fOut=fOut)

	print("Output: %s"%fOut)
	print("Size: (%d,%d,%d)"%(Ni,Nj,Nk))
	print("Inner Radius: %f"%Rin)
	print("Sunward Outer Radius: %f"%rOutS)
	print("Tail Outer Radius: %f"%rOutAS)
	print("Low-lat BC: %f"%(llBC))
	if (not doChimp):
		print("Ring params: \n%s"%(rParams[gid]))
	print("\nWriting to %s"%(fOut))

	if (doViz):
		gg.VizGrid(XX,YY,xxG,yyG,fOut=fOut,doGhost=doVizG)
	#gg.genRing(XX,YY,Nk=Nk,Tol=1.0,doVerb=True)
