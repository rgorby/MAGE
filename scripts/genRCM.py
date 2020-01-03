#!/usr/bin/env python
#Generates RCM config data

import argparse
import kaipy.rcm.rcminit as rcminit
from argparse import RawTextHelpFormatter
import numpy as np
import h5py

if __name__ == "__main__":

#Arg parsing
	fOut = "rcmconfig.h5"

	MainS = """Generates RCM configuration data

	"""
	
	parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
	parser.add_argument('-o',type=str,default=fOut,metavar="fOut",help="Output file name (default: %(default)s)")
	#Finalize parsing
	args = parser.parse_args()
	fOut = args.o
	print("Writing RCM configuration to %s"%(fOut))

	alamc,etac,ikflavc,fudgec = rcminit.LoadLAS1()
	dktab = rcminit.LoadDKTab()
	iflavin, alamin = rcminit.LoadEnchan()

	with h5py.File(fOut,'w') as hf:
		hf.create_dataset("alamc",data=alamc)
		hf.create_dataset("etac" ,data=etac)
		hf.create_dataset("ikflavc" ,data=ikflavc,dtype=np.int)
		hf.create_dataset("fudgec",data=fudgec)
		hf.create_dataset("dktable",data=dktab)
		hf.create_dataset("iflavin" ,data=iflavin,dtype=np.int)
		hf.create_dataset("alamin",data=alamin)
