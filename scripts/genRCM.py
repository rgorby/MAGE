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
        parser.add_argument('--p',metavar='True',choices=('True','False'),help="True if you want the no loss for the first channel")
	#Finalize parsing
	args = parser.parse_args()
	fOut = args.o
        p = args.p == 'True'
	print("Writing RCM configuration to %s"%(fOut))
        print("NOTE: This version only uses the file enchan.dat")

#	alamc,etac,ikflavc,fudgec = rcminit.LoadLAS1()
	dktab = rcminit.LoadDKTab()
 	ikflavc, alamc = rcminit.LoadEnchan()

        kdim = ikflavc.size
        print(' kdim = %d'%kdim)
        fudgec = np.zeros(kdim)
        for k in range(kdim):
         if(ikflavc[k] ==1):
          fudgec[k] = 0.33333
         else:
          fudgec[k] = 0.0

        if(p):
         print('Setting first channel to fudge =0')
         fudgec[0]=0.0

	with h5py.File(fOut,'w') as hf:
		hf.create_dataset("alamc",data=alamc)
		hf.create_dataset("ikflavc" ,data=ikflavc,dtype=np.int)
		hf.create_dataset("fudgec",data=fudgec)
		hf.create_dataset("dktable",data=dktab)
