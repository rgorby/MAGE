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
	""" If given single value, will return single datetime.datetime
		If given list, will return list of datetime.datetimes
	"""
	UT = Time(mjd,format='mjd').isot
	if type(UT) == str:
		return datetime.datetime.strptime(UT,isotfmt)
	else:
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

#transform from cartesian to spherical
def xyz2rtp(phi,theta,Ax,Ay,Az):
	Ar = Ax*np.cos(phi)*np.sin(theta)+Ay*np.sin(phi)*np.sin(theta)+Az*np.cos(theta)
	Ap = -Ax*np.sin(phi)+Ay*np.cos(phi)
	At = Ax*np.cos(phi)*np.cos(theta)+Ay*np.sin(phi)*np.cos(theta)-Az*np.sin(theta)
	return Ar,At,Ap
