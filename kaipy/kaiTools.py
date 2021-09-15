import numpy as np
import datetime
import os
import glob
import sys
import subprocess
from xml.dom import minidom

from astropy.time import Time

import h5py

isotfmt = '%Y-%m-%dT%H:%M:%S.%f'

def MJD2UT(mjd):
	""" Leaving in case new way breaks something. I did check and this method and new one give equivalent results but new one is a little faster
	astroT = Time(mjd,format='mjd').iso
	utall = []
	for ut in astroT:
		utall.append(datetime.datetime.strptime(ut,isotfmt))

	return utall
	"""
	UT = Time(mjd,format='mjd').isot
	return [datetime.datetime.strptime(UT[n],isotfmt) for n in range(len(UT))]

def getRunInfo(fdir,ftag):
	idStr = "_0000_0000_0000.gam.h5"
	isMPI = False
	Ri = 0
	Rj = 0
	Rk = 0
	fOld = os.path.join(fdir,ftag+'.h5')
	fNew = os.path.join(fdir,ftag+'.gam.h5')
	try:
		if (os.path.exists(fOld)):
			return fOld,isMPI,Ri,Rj,Rk
		if (os.path.exists(fNew)):
			return fNew,isMPI,Ri,Rj,Rk
		fIns = glob.glob(os.path.join(fdir,ftag)+'_*'+idStr)
		if (len(fIns) > 1):
			raise ValueError('Should not find more that one parallel file')
		if (len(fIns) == 0):
			raise ValueError('No MPI database found')
		else:
			isMPI = True
			fName = fIns[0]
			Ns = [int(s) for s in fName.split('_') if s.isdigit()]
			Ri = Ns[-5]
			Rj = Ns[-4]
			Rk = Ns[-3]
			return fName,isMPI,Ri,Rj,Rk
	except ValueError as ve:
			print(ve)
			sys.exit()

