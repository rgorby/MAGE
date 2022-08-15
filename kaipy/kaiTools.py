import numpy as np
import datetime
import os
import glob
import sys
import subprocess
from xml.dom import minidom
from astropy.time import Time
import h5py

import kaipy.kdefs as kd

isotfmt = '%Y-%m-%dT%H:%M:%S.%f'

to_center1D = lambda A: 0.5*(A[:-1]+A[1:])
to_center2D = lambda A: 0.25*(A[:-1,:-1]+A[1:,:-1]+A[:-1,1:]+A[1:,1:])
to_center3D = lambda A: 0.125*(A[:-1,:-1,:-1]+A[-1:,-1:,1:]+A[-1:,1:,-1:]+A[-1:,1:,1:]
				+A[1:,:-1,:-1]+A[1:,-1:,1:]+A[1:,1:,-1:]+A[1:,1:,1:])


def L_to_bVol(L, bsurf_nT=kd.EarthM0g*kd.G2nT):  # L shell [Rp] to V [Rp/nT]
	"""Calculates the flux tube volume [Rp/nT] from the given L shell [Rp]
		L : L shell [Rp]
		bsurf_nT : Surface field [nT] (Default: Earth)
	"""
	colat = np.arcsin(np.sqrt(1.0/L))

	cSum = 35*np.cos(colat) - 7*np.cos(3*colat) +(7./5.)*np.cos(5*colat) - (1./7.)*np.cos(7*colat)
	cSum /= 64.
	s8 = np.sin(colat)**8
	V = 2*cSum/s8/bsurf_nT
	return V  # [Rp/nT]

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

# Use the Burton 1975 Formula to compute Dst from solar wind parameters
def burtonDst(secs,n,vx,vy,vz,bx,by,bz):
	"""
	Given time in seconds, density in cm^-3, velocity in km/s and magnetic field in nT 
	this function will return Dst computed from the Burton 1975 (doi: 10.1029/ja080i031p04204)
	"""
	v = np.sqrt(vx**2+vy**2+vz**2)
	b = np.sqrt(bx**2+by**2+bz**2)
	ey = 1.0e-3*vx*bz
	a = 3.6e-5
	b = 0.2
	c = 20.0
	d = -1.5e-3
	dp = -1.2e-3
	fe = d*(np.maximum(0,ey-0.5))
	pdyn = 1.6726e-27*1e6*6.248e18*n*v**2
	sqrtpdyn = np.sqrt(pdyn)
	dst = 0*pdyn
	for i in np.arange(len(dst)-1)+1:
		dst[i] = dst[i-1] + (secs[i]-secs[i-1])*(fe[i]-a*(dst[i-1]-b*sqrtpdyn[i]+c))
	return dst

def newellkp(secs,n,vx,vy,vz,bx,by,bz):
	"""
	Given time in seconds, density in cm^-3, velocity in km/s and magnetic field in nT 
	this function will return Kp computed from the Newell 2008 (doi: 10.1029/2007ja012825)
	"""
	v = np.sqrt(vx**2+vy**2+vz**2)
	newellCoup = newellcoupling(vx,vy,vz,bx,by,bz)
	kp = np.clip(0.05+2.244e-4*newellCoup+2.844e-6*np.sqrt(n)*v**2,0,9.0)
	# Kp is  three hour index so compute a three hour rolling average using convolve with first and last value extended
	window = int(3*60*60/(secs[1]-secs[0]))
	halfwindow = int(window/2)
	begin = np.ones(halfwindow)
	end = np.ones(halfwindow)
	begin[:] = kp[0]
	end[:] = kp[-1]
	extendkp = np.concatenate((begin,kp,end))
	conextendkp = np.convolve(extendkp,np.ones(window)/window,mode='same')
	return conextendkp[halfwindow:-halfwindow]

def newellcoupling(vx,vy,vz,bx,by,bz):
	"""
	Given  density in cm^-3, velocity in km/s and magnetic field in nT 
	this function will return Universial Coupling Function computed from the Newell 2007 (doi: 10.1029/2006ja012015)
	"""
	v = np.sqrt(vx**2+vy**2+vz**2)
	b = np.sqrt(bx**2+by**2+bz**2)
	thetac = np.abs(np.arctan2(by,bz))
	newcoup = np.power(v,4.0/3.0)*np.power(b,2.0/3.0)*np.power(np.sin(thetac/2.0),8.0/3.0)
	return newcoup

#Read SymH from bcwind file
def GetSymH(fBC):
	with h5py.File(fBC,'r') as hf:
		mjdData = hf['MJD'][()]
		tData   = hf['T'][()]
		dstData = hf['symh'][()]
	
	utData = MJD2UT(mjdData)
	return utData,tData,dstData

def interpTSC(gridX, gridY, x, y, var):
	""" gridX/gridY: 3-element x & y grid vals (center value is closest grid point to desired interpolation point)
		x: dim1 of point of interest
		y: dim2 of point of interest
		var: 3x3 values of desired variable
	"""

	weights = interpTSCWeights(gridX, gridY, x, y)

	result = 0
	for i in range(3):
		for j in range(3):
			result += weights[i,j]*var[i,j]

	return result

def interpTSCWeights(gridX, gridY, x, y):
	""" gridX/gridY: 3-element x & y grid vals (center value is closest grid point to desired interpolation point)
		x: dim1 of point of interest
		y: dim2 of point of interest
	"""

	dx = np.abs(gridX[0]-gridX[1])
	dy = np.abs(gridY[0]-gridY[1])

	eta  = (x - gridX[1])/dx
	zeta = (y - gridY[1])/dy
	#print('eta : ' + str(eta))
	#print('zeta: ' + str(zeta))

	def weight1D(eta):
		return np.array([
				0.5*(0.5-eta)**2,
				0.75 - eta**2,
				0.5*(0.5+eta)**2
			])

	wX = weight1D(eta)
	wY = weight1D(zeta)

	w2D = np.zeros((3,3))
	for i in range(3):
		for j in range(3):
			w2D[i,j] = wX[i]*wY[j]

	return w2D
