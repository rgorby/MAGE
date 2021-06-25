#!/usr/bin/env python
#Takes MIX restart and up/down-scales it

import argparse
import kaipy.gamera.magsphereRescale as upscl
from argparse import RawTextHelpFormatter
import numpy as np
import h5py
import os
import kaipy.kaiH5 as kh5

if __name__ == "__main__":
	dIn = os.getcwd()


	inid  = "msphere"
	outid = "msphereX"

	nRes = "0"

	MainS = """Up/down-scales a RCM restart (kinda)

	inid/nres : Run ID string and restart number, i.e. input file = inid.MPISTUFF.Res.#nres.h5
	outid : Output Run ID

	"""	
	parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
	parser.add_argument('-i',metavar='inid',default=inid,help="Input Run ID string (default: %(default)s)")
	parser.add_argument('-n',type=int,metavar="nres",default=0,help="Restart number (default: %(default)s)")
	parser.add_argument('-o',type=str,metavar="outid",default=outid,help="Output run ID (default: %(default)s)")	
	parser.add_argument('--down',action='store_true',default=False,help='Downscale instead of upscale (default: %(default)s)')

	#Finalize parsing
	args = parser.parse_args()
	bStr = args.i
	nRes = args.n
	outid = args.o
	doUp = not args.down

#Start w/ mhd2imag restart
	fIn  = bStr  + ".mhd2imag.Res.%05d.h5"%(nRes)
	fOut = outid + ".mhd2imag.Res.%05d.h5"%(nRes)
	print("Reading from %s and writing to %s"%(fIn,fOut))

	vIDs = kh5.getRootVars(fIn)
	imNi,imNj = kh5.getDims(fIn,"Bmin") #Array size of coupler

	#Open input and output
	iH5 = h5py.File(fIn ,'r')
	oH5 = h5py.File(fOut,'w')

	#Start by scraping attributes
	for k in iH5.attrs.keys():
		aStr = str(k)
		#print(aStr)
		oH5.attrs.create(k,iH5.attrs[aStr])

	for vID in vIDs:
		Q = iH5[vID][:]
		
		if (doUp):
			Qr = upscl.upRCMCpl(Q,N=imNj)
		else:
			#Qr = upscl.downMIX(Q)
			print("Downscaling not implemented ...")
		oH5.create_dataset(vID,data=Qr)
		print("\t%s, Dims: %s => %s"%(vID,Q.shape,Qr.shape))

	#Now get 
	iH5.close()
	oH5.close()

#Now do RCM restart
	fIn  = bStr  + ".RCM.Res.%05d.h5"%(nRes)
	fOut = outid + ".RCM.Res.%05d.h5"%(nRes)
	print("\n\nReading from %s and writing to %s"%(fIn,fOut))

	vIDs = kh5.getRootVars(fIn)
	Ni,Nj,Nk = kh5.getDims(fIn,"rcmeeta")

	#Open input and output
	iH5 = h5py.File(fIn ,'r')
	oH5 = h5py.File(fOut,'w')

	#Start by scraping attributes
	for k in iH5.attrs.keys():
		aStr = str(k)
		#print(aStr)
		oH5.attrs.create(k,iH5.attrs[aStr])

	for vID in vIDs:
		Q = iH5[vID][:]
		
		if (doUp):
			Qr = upscl.upRCM(Q,Ni=Ni,Nj=Nj,Nk=Nk)
		else:
			#Qr = upscl.downMIX(Q)
			print("Downscaling not implemented ...")
		print("\t%s, Dims: %s => %s"%(vID,Q.shape,Qr.shape))
		oH5.create_dataset(vID,data=Qr)
		
	#Now finish
	iH5.close()
	oH5.close()

	Nri = Ni
	Nrj = (Nj-2)*2+2
	print("\n\nFinished rescaling, new RCM config must be:")
	print("\tRCMSIZEI = %d"%(Nri))
	print("\tRCMSIZEJ = %d"%(Nrj))
	print("\tRCMSIZEK = %d"%(Nk))
